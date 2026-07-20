#include "kmat/pca.hpp"

#include "kmat/io.hpp"
#include "kmat/table.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <random>

namespace kmat {

Error run_pca(const PcaOptions& opts, PcaResult& result) {
  std::vector<std::string> accession_paths;
  if (auto err = read_list_file(opts.accession_list_path, accession_paths); !err.ok()) {
    return err;
  }
  if (accession_paths.empty()) {
    return Error::invalid_argument("accession list is empty");
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
    if (matrix.header.num_accessions != accession_paths.size()) {
      return Error::invalid_argument("matrix accession count mismatch");
    }
  } else {
    return Error::invalid_argument("matrix path or matrix list required");
  }

  if (matrix.patterns.empty()) {
    return Error::invalid_argument("matrix has no patterns");
  }

  std::vector<std::size_t> sample_indices(matrix.patterns.size());
  for (std::size_t i = 0; i < sample_indices.size(); ++i) {
    sample_indices[i] = i;
  }

  if (opts.max_samples > 0 && opts.max_samples < sample_indices.size()) {
    std::mt19937 rng(opts.seed);
    std::shuffle(sample_indices.begin(), sample_indices.end(), rng);
    sample_indices.resize(opts.max_samples);
    std::sort(sample_indices.begin(), sample_indices.end());
  }

  const std::size_t n_acc = accession_paths.size();
  const std::size_t n_samples = sample_indices.size();
  const std::size_t npc = std::min(opts.num_pcs, std::min(n_acc, n_samples));

  if (npc == 0) {
    return Error::invalid_argument("cannot compute zero principal components");
  }

  Eigen::MatrixXd X(static_cast<Eigen::Index>(n_acc), static_cast<Eigen::Index>(n_samples));
  for (std::size_t s = 0; s < n_samples; ++s) {
    const auto& words = matrix.patterns[sample_indices[s]];
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
