#pragma once

#include "kmat/error.hpp"
#include "kmat/presence.hpp"

#include <cstddef>
#include <string>

namespace kmat {

struct CountOptions {
  std::string input_path;
  std::string output_path;
  std::size_t kmer_size{31};
  std::uint32_t min_count{1};  // KMC-style -ci
};

/// Count k-mers from FASTQ/FASTA (.gz ok); keep codes with count >= min_count; write `.kset`.
Error count_kmers_to_presence_set(const CountOptions& opts);

struct ImportKmersOptions {
  std::string input_path;   // one DNA k-mer per line (or hex uint64 with 0x prefix)
  std::string output_path;
  std::size_t kmer_size{31};
};

/// Migration helper without linking KMC: import a text k-mer list into a `.kset`.
Error import_kmers_text_to_presence_set(const ImportKmersOptions& opts);

}  // namespace kmat
