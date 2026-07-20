#pragma once

#include "kmat/error.hpp"

#include <string>
#include <vector>

namespace kmat {

struct FastaRecord {
  std::string id;
  std::string sequence;
};

/// Read all sequences from a FASTA file (concatenates multi-record files).
Error read_fasta_sequences(const std::string& path, std::vector<FastaRecord>& records);

/// Uppercase A/C/G/T/N stripped; other chars skipped.
Error normalize_dna(std::string& sequence);

}  // namespace kmat
