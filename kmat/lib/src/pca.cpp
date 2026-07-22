#include "kmat/pca.hpp"

#include "kmat/io.hpp"
#include "kmat/log.hpp"
#include "kmat/table.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <fstream>
#include <random>
#include <unordered_set>

namespace kmat {

namespace {

constexpr std::size_t kDefaultMaxSamples = 100000;
constexpr std::size_t kAutoCapPatterns = 500000;

std::vector<std::size_t> sample_unique_indices(std::size_t n, std::size_t k, unsigned seed) {
  std::vector<std::size_t> out;
  if (k == 0 || n == 0) {
    return out;
  }
  if (k >= n) {
    out.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      out[i] = i;
    }
    return out;
  }

  std::mt19937 rng(seed);
  // For moderate k, draw until unique; for large k relative to n, shuffle a partial range.
  if (k * 4 < n) {
    std::unordered_set<std::size_t> seen;
    seen.reserve(k * 2);
    std::uniform_int_distribution<std::size_t> dist(0, n - 1);
    while (seen.size() < k) {
      seen.insert(dist(rng));
    }
    out.assign(seen.begin(), seen.end());
    std::sort(out.begin(), out.end());
    return out;
  }

  out.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = i;
  }
  std::shuffle(out.begin(), out.end(), rng);
  out.resize(k);
  std::sort(out.begin(), out.end());
  return out;
}

/// Load only the pattern store from a v2 matrix (skip the huge k-mer map).
Error read_v2_patterns_only(const std::string& path, MatrixHeader& header,
                            std::vector<std::vector<std::uint64_t>>& patterns) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return Error::io_error("failed to open matrix: " + path);
  }
  in.read(reinterpret_cast<char*>(&header), sizeof(header));
  if (!in) {
    return Error::io_error("failed reading matrix header");
  }
  if (header.magic[0] != 'K' || header.magic[1] != 'M' || header.magic[2] != 'A' ||
      header.magic[3] != 'T') {
    return Error::invalid_argument("invalid matrix magic");
  }
  if (header.version != 2) {
    return Error::invalid_argument("pop streaming path requires v2 matrix");
  }

  const std::uint64_t num_patterns = header.reserved;
  patterns.resize(static_cast<std::size_t>(num_patterns));
  for (std::size_t p = 0; p < patterns.size(); ++p) {
    patterns[p].resize(header.num_stripes);
    in.read(reinterpret_cast<char*>(patterns[p].data()),
            static_cast<std::streamsize>(header.num_stripes * sizeof(std::uint64_t)));
    if (!in) {
      return Error::io_error("unexpected EOF reading pattern store");
    }
  }
  return Error::success();
}

/// Reservoir-sample patterns while reading so peak RAM stays O(max_samples), not O(all patterns).
Error reservoir_sample_v2_patterns(const std::string& path, std::size_t max_samples, unsigned seed,
                                   MatrixHeader& header,
                                   std::vector<std::vector<std::uint64_t>>& sampled) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return Error::io_error("failed to open matrix: " + path);
  }
  in.read(reinterpret_cast<char*>(&header), sizeof(header));
  if (!in) {
    return Error::io_error("failed reading matrix header");
  }
  if (header.magic[0] != 'K' || header.magic[1] != 'M' || header.magic[2] != 'A' ||
      header.magic[3] != 'T') {
    return Error::invalid_argument("invalid matrix magic");
  }
  if (header.version != 2) {
    return Error::invalid_argument("pop reservoir path requires v2 matrix");
  }

  const std::uint64_t num_patterns = header.reserved;
  if (num_patterns == 0) {
    return Error::invalid_argument("matrix has no patterns");
  }

  std::mt19937 rng(seed);
  sampled.clear();
  sampled.reserve(std::min<std::size_t>(max_samples, static_cast<std::size_t>(num_patterns)));

  std::vector<std::uint64_t> words(header.num_stripes);
  for (std::uint64_t i = 0; i < num_patterns; ++i) {
    in.read(reinterpret_cast<char*>(words.data()),
            static_cast<std::streamsize>(header.num_stripes * sizeof(std::uint64_t)));
    if (!in) {
      return Error::io_error("unexpected EOF reading pattern store");
    }
    if (sampled.size() < max_samples) {
      sampled.push_back(words);
    } else {
      std::uniform_int_distribution<std::uint64_t> dist(0, i);
      const std::uint64_t j = dist(rng);
      if (j < max_samples) {
        sampled[static_cast<std::size_t>(j)] = words;
      }
    }
  }
  return Error::success();
}

}  // namespace

Error run_pca(const PcaOptions& opts, PcaResult& result) {
  std::vector<std::string> accession_paths;
  if (auto err = read_list_file(opts.accession_list_path, accession_paths); !err.ok()) {
    return err;
  }
  if (accession_paths.empty()) {
    return Error::invalid_argument("accession list is empty");
  }

  const std::size_t n_acc = accession_paths.size();
  std::vector<std::vector<std::uint64_t>> patterns;
  MatrixHeader header{};

  if (!opts.matrix_list_path.empty()) {
    PaMatrix matrix;
    if (auto err = load_matrix_from_list(opts.matrix_list_path, accession_paths.size(), matrix);
        !err.ok()) {
      return err;
    }
    header = matrix.header;
    patterns = std::move(matrix.patterns);
  } else if (!opts.matrix_path.empty()) {
    // Peek header to choose streaming vs full-pattern load.
    {
      std::ifstream in(opts.matrix_path, std::ios::binary);
      if (!in) {
        return Error::io_error("failed to open matrix: " + opts.matrix_path);
      }
      in.read(reinterpret_cast<char*>(&header), sizeof(header));
      if (!in) {
        return Error::io_error("failed reading matrix header");
      }
    }

    std::size_t max_samples = opts.max_samples;
    if (max_samples == 0) {
      if (header.reserved > kAutoCapPatterns) {
        max_samples = kDefaultMaxSamples;
        log_warn("pop: auto-capping pattern samples to " + std::to_string(max_samples) +
                 " (use --max-samples to override; 0 no longer means all on large panels)");
      } else {
        max_samples = static_cast<std::size_t>(header.reserved);
      }
    }

    if (header.version == 2 && max_samples < header.reserved) {
      log_info("pop: reservoir-sampling " + std::to_string(max_samples) + " / " +
               std::to_string(header.reserved) + " patterns (skipping k-mer map)");
      if (auto err = reservoir_sample_v2_patterns(opts.matrix_path, max_samples, opts.seed, header,
                                                  patterns);
          !err.ok()) {
        return err;
      }
    } else if (header.version == 2) {
      log_info("pop: loading all " + std::to_string(header.reserved) +
               " patterns (skipping k-mer map)");
      if (auto err = read_v2_patterns_only(opts.matrix_path, header, patterns); !err.ok()) {
        return err;
      }
    } else {
      PaMatrix matrix;
      if (auto err = read_matrix(opts.matrix_path, matrix); !err.ok()) {
        return err;
      }
      header = matrix.header;
      patterns = std::move(matrix.patterns);
    }

    if (header.num_accessions != accession_paths.size()) {
      return Error::invalid_argument("matrix accession count mismatch");
    }
  } else {
    return Error::invalid_argument("matrix path or matrix list required");
  }

  if (patterns.empty()) {
    return Error::invalid_argument("matrix has no patterns");
  }

  std::size_t max_samples = opts.max_samples;
  if (max_samples == 0) {
    if (patterns.size() > kAutoCapPatterns) {
      max_samples = kDefaultMaxSamples;
    } else {
      max_samples = patterns.size();
    }
  }

  std::vector<std::size_t> sample_indices;
  if (max_samples >= patterns.size()) {
    sample_indices.resize(patterns.size());
    for (std::size_t i = 0; i < patterns.size(); ++i) {
      sample_indices[i] = i;
    }
  } else if (patterns.size() == max_samples) {
    // Already reservoir-sampled into `patterns`.
    sample_indices.resize(patterns.size());
    for (std::size_t i = 0; i < patterns.size(); ++i) {
      sample_indices[i] = i;
    }
  } else {
    sample_indices = sample_unique_indices(patterns.size(), max_samples, opts.seed);
  }

  const std::size_t n_samples = sample_indices.size();
  const std::size_t npc_cap = std::min(n_acc, n_samples);
  const std::size_t npc =
      (opts.num_pcs == 0) ? npc_cap : std::min(opts.num_pcs, npc_cap);

  if (npc == 0) {
    return Error::invalid_argument("cannot compute zero principal components");
  }

  log_info("pop: PCA accessions=" + std::to_string(n_acc) + " samples=" +
           std::to_string(n_samples) + " npc=" + std::to_string(npc));

  Eigen::MatrixXd X(static_cast<Eigen::Index>(n_acc), static_cast<Eigen::Index>(n_samples));
  for (std::size_t s = 0; s < n_samples; ++s) {
    const auto& words = patterns[sample_indices[s]];
    for (std::size_t a = 0; a < n_acc; ++a) {
      X(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(s)) =
          get_presence_bit(words, a) ? 1.0 : 0.0;
    }
  }

  Eigen::MatrixXd centered = X.rowwise() - X.colwise().mean();
  Eigen::MatrixXd cov =
      (centered * centered.transpose()) / static_cast<double>(std::max<std::size_t>(1, n_samples - 1));

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cov);
  if (solver.info() != Eigen::Success) {
    return Error::invalid_argument("PCA eigen decomposition failed");
  }

  const Eigen::MatrixXd evecs = solver.eigenvectors();
  const std::size_t start_col = (n_acc >= npc) ? (n_acc - npc) : 0;

  result.accessions.clear();
  result.accessions.reserve(n_acc);
  for (const std::string& path : accession_paths) {
    result.accessions.push_back(accession_id_from_path(path));
  }

  result.pcs.assign(n_acc, std::vector<double>(npc, 0.0));
  for (std::size_t a = 0; a < n_acc; ++a) {
    for (std::size_t pc = 0; pc < npc; ++pc) {
      result.pcs[a][pc] = evecs(static_cast<Eigen::Index>(a),
                                start_col + static_cast<Eigen::Index>(pc));
    }
  }

  return Error::success();
}

Error write_pca_tsv(const std::string& path, const PcaResult& result) {
  return write_pop_tsv(path, result.accessions, result.pcs);
}

}  // namespace kmat
