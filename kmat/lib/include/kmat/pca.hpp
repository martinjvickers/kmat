#pragma once

#include "kmat/error.hpp"
#include "kmat/matrix.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace kmat {

struct PcaOptions {
  std::string matrix_path;
  std::string matrix_list_path;
  std::string accession_list_path;
  std::size_t num_pcs{2};
  std::size_t max_samples{0};  // 0 = all rows
  unsigned seed{42};
};

struct PcaResult {
  std::vector<std::string> accessions;
  std::vector<std::vector<double>> pcs;
};

Error run_pca(const PcaOptions& opts, PcaResult& result);

Error write_pca_tsv(const std::string& path, const PcaResult& result);

}  // namespace kmat
