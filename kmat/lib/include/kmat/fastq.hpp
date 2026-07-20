#pragma once

#include "kmat/error.hpp"
#include "kmat/fasta.hpp"

#include <string>
#include <vector>

namespace kmat {

struct FastqRecord {
  std::string id;
  std::string sequence;
  std::string quality;
};

/// Read all records from a FASTQ file (plain or .gz).
Error read_fastq_sequences(const std::string& path, std::vector<FastqRecord>& records);

/// Load DNA sequences from FASTA/FASTQ, optionally gzip-compressed (by extension).
Error read_sequence_file(const std::string& path, std::vector<FastaRecord>& records);

bool path_looks_gzip(const std::string& path);
bool path_looks_fastq(const std::string& path);
bool path_looks_fasta(const std::string& path);

}  // namespace kmat
