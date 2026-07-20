#pragma once

#include "kmat/error.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace kmat {

struct GeneSearchOptions {
  std::string matrix_path;
  std::string matrix_list_path;
  std::string accession_list_path;
  std::string gene_fasta_path;
  std::size_t kmer_size{31};
};

struct GeneHit {
  std::uint64_t kmer_code{0};
  std::string kmer;
  std::string pa_bits;
};

struct GeneSearchResult {
  std::vector<GeneHit> hits;
};

Error run_gene_search(const GeneSearchOptions& opts, GeneSearchResult& result);

Error write_gene_hits(std::ostream& out, const GeneSearchResult& result);

}  // namespace kmat
