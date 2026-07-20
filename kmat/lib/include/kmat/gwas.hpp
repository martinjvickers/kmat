#pragma once

#include "kmat/error.hpp"
#include "kmat/matrix.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace kmat {

struct GwasOptions {
  std::string matrix_path;
  std::string matrix_list_path;
  std::string accession_list_path;
  std::string phenotype_path;
  std::string pop_path;
  std::size_t kmer_size{31};
  std::size_t top_n{1000};
  bool print_all{false};
  /// Fill `pa_bits` on hits (expensive); default off.
  bool include_pa_bits{false};
  /// 0 = use process runtime profile thread count.
  std::size_t num_threads{0};
};

struct GwasHit {
  std::uint64_t kmer_code{0};
  std::vector<double> phenotype_sums;
  std::size_t num_carriers{0};
  std::vector<double> p_values;
  std::string pa_bits;
};

struct GwasResult {
  std::vector<std::string> phenotype_names;
  std::vector<GwasHit> hits;
};

Error run_gwas(const GwasOptions& opts, GwasResult& result);

Error write_gwas_tsv(std::ostream& out, std::size_t kmer_size, const GwasResult& result,
                     bool include_pa_bits);

}  // namespace kmat
