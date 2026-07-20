#pragma once

#include "kmat/error.hpp"

#include <cstdint>
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

bool path_looks_presence_set(const std::string& path);

}  // namespace kmat
