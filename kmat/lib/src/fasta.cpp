#include "kmat/fasta.hpp"

#include <cctype>
#include <fstream>

namespace kmat {

Error normalize_dna(std::string& sequence) {
  std::string cleaned;
  cleaned.reserve(sequence.size());
  for (char c : sequence) {
    const char u = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (u == 'A' || u == 'C' || u == 'G' || u == 'T' || u == 'N') {
      cleaned.push_back(u);
    }
  }
  sequence.swap(cleaned);
  return Error::success();
}

Error read_fasta_sequences(const std::string& path, std::vector<FastaRecord>& records) {
  std::ifstream in(path);
  if (!in) {
    return Error::io_error("failed to open FASTA: " + path);
  }

  records.clear();
  FastaRecord current;
  bool have_record = false;
  std::string line;

  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty()) {
      continue;
    }
    if (line.front() == '>') {
      if (have_record) {
        if (auto err = normalize_dna(current.sequence); !err.ok()) {
          return err;
        }
        records.push_back(std::move(current));
        current = FastaRecord{};
      }
      current.id = line.substr(1);
      have_record = true;
      continue;
    }
    if (!have_record) {
      return Error::invalid_argument("FASTA sequence before header in " + path);
    }
    current.sequence += line;
  }

  if (have_record) {
    if (auto err = normalize_dna(current.sequence); !err.ok()) {
      return err;
    }
    records.push_back(std::move(current));
  }

  if (records.empty()) {
    return Error::invalid_argument("no sequences found in FASTA: " + path);
  }
  return Error::success();
}

}  // namespace kmat
