#pragma once

#include "kmat/error.hpp"

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace kmat {

/// On-disk presence set (`.kset`): filtered k-mers for one accession.
struct PresenceHeader {
  char magic[4]{'K', 'S', 'E', 'T'};
  std::uint32_t version{1};
  std::uint32_t kmer_size{0};
  std::uint32_t min_count{1};
  std::uint32_t reserved{0};
  std::uint64_t num_kmers{0};
};

static_assert(sizeof(PresenceHeader) == 32, "PresenceHeader must be 32 bytes");

struct PresenceSet {
  PresenceHeader header{};
  std::vector<std::uint64_t> kmers;  // sorted unique codes
};

Error write_presence_set(const std::string& path, const PresenceSet& set);

Error read_presence_set(const std::string& path, PresenceSet& set);

/// Stream sorted k-mer codes from a `.kset` without loading the full vector.
class PresenceSetCursor {
 public:
  Error open(const std::string& path);
  void close();

  [[nodiscard]] bool ok() const { return ok_; }
  [[nodiscard]] bool has_value() const { return has_value_; }
  [[nodiscard]] std::uint64_t value() const { return value_; }
  [[nodiscard]] const PresenceHeader& header() const { return header_; }
  [[nodiscard]] const std::string& path() const { return path_; }

  /// Advance to the next k-mer. Clears has_value at EOF.
  Error advance();

 private:
  std::string path_;
  std::ifstream in_;
  PresenceHeader header_{};
  std::uint64_t remaining_{0};
  std::uint64_t value_{0};
  bool has_value_{false};
  bool ok_{false};
};

bool path_looks_presence_set(const std::string& path);

}  // namespace kmat
