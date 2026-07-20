#include "kmat/count.hpp"

#include "kmat/fastq.hpp"
#include "kmat/kmer.hpp"
#include "kmat/sequence.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace kmat {

Error count_kmers_to_presence_set(const CountOptions& opts) {
  if (opts.input_path.empty() || opts.output_path.empty()) {
    return Error::invalid_argument("count requires input and output paths");
  }
  if (opts.kmer_size == 0 || opts.kmer_size > 32) {
    return Error::invalid_argument("k-mer size must be in 1..32");
  }
  if (opts.min_count == 0) {
    return Error::invalid_argument("min_count (--ci) must be >= 1");
  }

  std::vector<FastaRecord> records;
  if (auto err = read_sequence_file(opts.input_path, records); !err.ok()) {
    return err;
  }

  std::unordered_map<std::uint64_t, std::uint32_t> counts;
  counts.reserve(1 << 16);
  for (const FastaRecord& rec : records) {
    if (auto err = for_each_encoded_kmer(rec.sequence, opts.kmer_size, [&](std::uint64_t code) {
          auto& c = counts[code];
          if (c < std::numeric_limits<std::uint32_t>::max()) {
            ++c;
          }
        });
        !err.ok()) {
      return err;
    }
  }

  PresenceSet set;
  set.header.kmer_size = static_cast<std::uint32_t>(opts.kmer_size);
  set.header.min_count = opts.min_count;
  set.kmers.reserve(counts.size());
  for (const auto& [code, c] : counts) {
    if (c >= opts.min_count) {
      set.kmers.push_back(code);
    }
  }
  std::sort(set.kmers.begin(), set.kmers.end());
  set.header.num_kmers = set.kmers.size();
  return write_presence_set(opts.output_path, set);
}

Error import_kmers_text_to_presence_set(const ImportKmersOptions& opts) {
  if (opts.input_path.empty() || opts.output_path.empty()) {
    return Error::invalid_argument("import-kmers requires input and output paths");
  }
  if (opts.kmer_size == 0 || opts.kmer_size > 32) {
    return Error::invalid_argument("k-mer size must be in 1..32");
  }

  std::ifstream in(opts.input_path);
  if (!in) {
    return Error::io_error("failed to open k-mer list: " + opts.input_path);
  }

  std::unordered_set<std::uint64_t> uniq;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty() || line.front() == '#') {
      continue;
    }
    while (!line.empty() && std::isspace(static_cast<unsigned char>(line.front()))) {
      line.erase(line.begin());
    }
    while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back()))) {
      line.pop_back();
    }
    if (line.empty()) {
      continue;
    }

    std::uint64_t code = 0;
    if (line.size() >= 2 && line[0] == '0' && (line[1] == 'x' || line[1] == 'X')) {
      try {
        code = std::stoull(line, nullptr, 16);
      } catch (...) {
        return Error::invalid_argument("invalid hex k-mer: " + line);
      }
    } else {
      if (line.size() != opts.kmer_size) {
        return Error::invalid_argument("k-mer length mismatch for: " + line);
      }
      if (auto err = encode_kmer(line, code); !err.ok()) {
        return err;
      }
    }
    uniq.insert(code);
  }

  PresenceSet set;
  set.header.kmer_size = static_cast<std::uint32_t>(opts.kmer_size);
  set.header.min_count = 1;
  set.kmers.assign(uniq.begin(), uniq.end());
  std::sort(set.kmers.begin(), set.kmers.end());
  set.header.num_kmers = set.kmers.size();
  return write_presence_set(opts.output_path, set);
}

}  // namespace kmat
