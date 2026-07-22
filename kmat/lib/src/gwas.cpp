#include "kmat/gwas.hpp"

#include "kmat/io.hpp"
#include "kmat/kmer.hpp"
#include "kmat/runtime.hpp"
#include "kmat/stats.hpp"
#include "kmat/table.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <iomanip>
#include <map>
#include <optional>

namespace kmat {

namespace {

Error accession_ids_from_paths(const std::vector<std::string>& paths,
                               std::vector<std::string>& ids) {
  ids.clear();
  ids.reserve(paths.size());
  for (const std::string& path : paths) {
    ids.push_back(accession_id_from_path(path));
  }
  return Error::success();
}

Error load_pop_table(const std::string& path, std::map<std::string, std::vector<double>>& pop) {
  std::vector<std::string> header;
  std::vector<TableRow> rows;
  if (auto err = read_numeric_table(path, header, rows); !err.ok()) {
    return err;
  }
  pop.clear();
  for (const TableRow& row : rows) {
    if (row.values.empty()) {
      return Error::invalid_argument("population structure row has no PC columns: " + row.key);
    }
    pop[row.key] = row.values;
  }
  return Error::success();
}

}  // namespace

Error run_gwas(const GwasOptions& opts, GwasResult& result) {
  std::vector<std::string> accession_paths;
  if (auto err = read_list_file(opts.accession_list_path, accession_paths); !err.ok()) {
    return err;
  }

  PaMatrix matrix;
  if (!opts.matrix_list_path.empty()) {
    if (auto err = load_matrix_from_list(opts.matrix_list_path, accession_paths.size(), matrix);
        !err.ok()) {
      return err;
    }
  } else if (!opts.matrix_path.empty()) {
    if (auto err = read_matrix(opts.matrix_path, matrix); !err.ok()) {
      return err;
    }
  } else {
    return Error::invalid_argument("matrix path or matrix list required");
  }

  if (matrix.header.kmer_size != opts.kmer_size) {
    return Error::invalid_argument("k-mer size mismatch between CLI and matrix");
  }

  std::vector<std::string> accession_ids;
  if (auto err = accession_ids_from_paths(accession_paths, accession_ids); !err.ok()) {
    return err;
  }

  std::vector<std::string> pheno_header;
  std::vector<TableRow> pheno_rows;
  if (auto err = read_numeric_table(opts.phenotype_path, pheno_header, pheno_rows); !err.ok()) {
    return err;
  }

  std::map<std::string, std::vector<double>> pop_map;
  if (auto err = load_pop_table(opts.pop_path, pop_map); !err.ok()) {
    return err;
  }

  std::map<std::string, std::size_t> accession_index;
  for (std::size_t i = 0; i < accession_ids.size(); ++i) {
    accession_index[accession_ids[i]] = i;
  }

  std::vector<int> pheno_to_accession;
  pheno_to_accession.reserve(pheno_rows.size());
  for (const TableRow& row : pheno_rows) {
    const auto it = accession_index.find(row.key);
    if (it == accession_index.end()) {
      return Error::invalid_argument("phenotype accession not in accession list: " + row.key);
    }
    pheno_to_accession.push_back(static_cast<int>(it->second));
  }

  const int N = static_cast<int>(pheno_to_accession.size());
  const int P = static_cast<int>(pheno_rows.front().values.size());
  if (opts.num_pcs == 0) {
    return Error::invalid_argument("gwas --npc must be >= 1");
  }
  const int npc = static_cast<int>(opts.num_pcs);
  const int n_cov = 1 + npc;
  const int df = N - n_cov;
  if (df <= 0) {
    return Error::invalid_argument("not enough phenotyped accessions for GWAS (need N > " +
                                   std::to_string(n_cov) + " for --npc " + std::to_string(npc) +
                                   ")");
  }

  result.phenotype_names.assign(pheno_header.begin() + 1, pheno_header.end());

  Eigen::MatrixXd Y(N, P);
  Eigen::MatrixXd C(N, n_cov);
  std::vector<int> word_idx(static_cast<std::size_t>(N));
  std::vector<std::uint64_t> bit_mask(static_cast<std::size_t>(N));

  for (int r = 0; r < N; ++r) {
    const int acc_i = pheno_to_accession[static_cast<std::size_t>(r)];
    const std::string& acc_id = accession_ids[static_cast<std::size_t>(acc_i)];

    for (int v = 0; v < P; ++v) {
      Y(r, v) = pheno_rows[static_cast<std::size_t>(r)].values[static_cast<std::size_t>(v)];
    }

    const auto pop_it = pop_map.find(acc_id);
    if (pop_it == pop_map.end()) {
      return Error::invalid_argument("population structure missing accession: " + acc_id);
    }
    if (pop_it->second.size() < opts.num_pcs) {
      return Error::invalid_argument(
          "population structure for " + acc_id + " has " + std::to_string(pop_it->second.size()) +
          " PC columns but --npc=" + std::to_string(opts.num_pcs));
    }
    C(r, 0) = 1.0;
    for (int pc = 0; pc < npc; ++pc) {
      C(r, 1 + pc) = pop_it->second[static_cast<std::size_t>(pc)];
    }

    word_idx[static_cast<std::size_t>(r)] = acc_i >> 6;
    bit_mask[static_cast<std::size_t>(r)] = (1ULL << (acc_i & 63));
  }

  // Residualize phenotypes once (shared across all patterns / threads).
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(C);
  Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(N, n_cov);
  Eigen::MatrixXd Y_res = Y - Q * (Q.transpose() * Y);
  Eigen::RowVectorXd Syy_res = Y_res.colwise().squaredNorm();

  struct ScoredPattern {
    double min_p;
    GwasHit template_hit;
    std::uint32_t pattern_id{0};
  };

  const std::size_t n_patterns = matrix.patterns.size();
  std::vector<std::optional<ScoredPattern>> slot(n_patterns);

  const std::size_t threads =
      opts.num_threads > 0 ? opts.num_threads : effective_threads(runtime_config());

  parallel_for(0, n_patterns, threads, [&](std::size_t pid_sz) {
    const std::uint32_t pattern_id = static_cast<std::uint32_t>(pid_sz);
    const auto& words = matrix.patterns[pattern_id];

    Eigen::VectorXd z(N);
    Eigen::VectorXd z_res(N);
    Eigen::VectorXd vec_eig(P);
    Eigen::RowVectorXd g(P);
    Eigen::VectorXd z_proj(n_cov);

    int num_bits_set = 0;
    for (int r = 0; r < N; ++r) {
      const bool present =
          (words[static_cast<std::size_t>(word_idx[static_cast<std::size_t>(r)])] &
           bit_mask[static_cast<std::size_t>(r)]) != 0ULL;
      z(r) = present ? 1.0 : 0.0;
      num_bits_set += present ? 1 : 0;
    }
    if (num_bits_set == 0) {
      return;
    }

    vec_eig.noalias() = Y.transpose() * z;
    z_proj.noalias() = Q.transpose() * z;
    z_res.noalias() = z - Q * z_proj;
    const double szz = z_res.squaredNorm();
    if (szz == 0.0) {
      return;
    }

    g.noalias() = z_res.transpose() * Y_res;
    const Eigen::RowVectorXd beta = g.array() / szz;
    const Eigen::RowVectorXd RSS = Syy_res.array() - beta.array() * g.array();
    const Eigen::RowVectorXd sigma2 = RSS.array() / static_cast<double>(df);
    const Eigen::RowVectorXd se = (sigma2.array() / szz).sqrt();
    const Eigen::RowVectorXd trow = beta.array() / se.array();

    GwasHit hit;
    hit.kmer_code = 0;
    hit.num_carriers = static_cast<std::size_t>(num_bits_set);
    hit.phenotype_sums.resize(static_cast<std::size_t>(P));
    hit.p_values.resize(static_cast<std::size_t>(P));
    for (int v = 0; v < P; ++v) {
      hit.phenotype_sums[static_cast<std::size_t>(v)] = vec_eig(v);
      hit.p_values[static_cast<std::size_t>(v)] = student_t_pvalue_two_sided(trow(v), df);
    }
    if (opts.include_pa_bits) {
      hit.pa_bits = presence_bitstring(words, accession_paths.size());
    }

    const double min_p = *std::min_element(hit.p_values.begin(), hit.p_values.end());
    slot[pid_sz] = ScoredPattern{min_p, std::move(hit), pattern_id};
  });

  std::vector<ScoredPattern> scored_patterns;
  scored_patterns.reserve(n_patterns);
  for (auto& s : slot) {
    if (s) {
      scored_patterns.push_back(std::move(*s));
    }
  }

  std::sort(scored_patterns.begin(), scored_patterns.end(),
            [](const ScoredPattern& a, const ScoredPattern& b) { return a.min_p < b.min_p; });

  if (!opts.print_all && scored_patterns.size() > opts.top_n) {
    scored_patterns.resize(opts.top_n);
  }

  const auto kmer_index = pattern_kmer_index(matrix);

  result.hits.clear();
  for (const ScoredPattern& scored : scored_patterns) {
    const auto& codes = kmer_index[scored.pattern_id];
    for (std::uint64_t code : codes) {
      GwasHit hit = scored.template_hit;
      hit.kmer_code = code;
      result.hits.push_back(std::move(hit));
    }
  }
  std::sort(result.hits.begin(), result.hits.end(), [](const GwasHit& a, const GwasHit& b) {
    const double min_a = *std::min_element(a.p_values.begin(), a.p_values.end());
    const double min_b = *std::min_element(b.p_values.begin(), b.p_values.end());
    if (min_a != min_b) {
      return min_a < min_b;
    }
    return a.kmer_code < b.kmer_code;
  });

  return Error::success();
}

Error write_gwas_tsv(std::ostream& out, std::size_t kmer_size, const GwasResult& result,
                     bool include_pa_bits) {
  out << "kmer";
  for (const std::string& name : result.phenotype_names) {
    out << '\t' << name << "_sum";
  }
  out << "\tnum_carriers";
  for (const std::string& name : result.phenotype_names) {
    out << '\t' << name << "_p";
  }
  if (include_pa_bits) {
    out << "\tpa_bits";
  }
  out << '\n';

  out << std::scientific << std::setprecision(6);
  for (const GwasHit& hit : result.hits) {
    std::string kmer;
    if (auto err = decode_kmer(hit.kmer_code, kmer_size, kmer); !err.ok()) {
      return err;
    }
    out << kmer;
    for (double v : hit.phenotype_sums) {
      out << '\t' << v;
    }
    out << '\t' << hit.num_carriers;
    for (double p : hit.p_values) {
      out << '\t' << p;
    }
    if (include_pa_bits) {
      out << '\t' << hit.pa_bits;
    }
    out << '\n';
  }
  out << std::defaultfloat;
  return Error::success();
}

}  // namespace kmat
