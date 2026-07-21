#pragma once

#include "kmat/error.hpp"

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace kmat {

/// On-disk master universe (`.kuniv`): sorted unique k-mer codes for the panel.
struct UniverseHeader {
  char magic[4]{'K', 'U', 'N', 'I'};
  std::uint32_t version{1};
  std::uint32_t kmer_size{0};
  std::uint32_t reserved{0};
  std::uint64_t num_kmers{0};
};

static_assert(sizeof(UniverseHeader) == 24, "UniverseHeader must be 24 bytes");

struct UniverseSet {
  UniverseHeader header{};
  std::vector<std::uint64_t> kmers;
};

Error write_universe(const std::string& path, const UniverseSet& set);

Error read_universe_header(const std::string& path, UniverseHeader& header);

/// Stream sorted unique codes from a `.kuniv` without loading the full vector.
class UniverseCursor {
 public:
  Error open(const std::string& path);
  void close();

  [[nodiscard]] bool has_value() const { return has_value_; }
  [[nodiscard]] std::uint64_t value() const { return value_; }
  [[nodiscard]] const UniverseHeader& header() const { return header_; }

  Error advance();

 private:
  std::ifstream in_;
  UniverseHeader header_{};
  std::uint64_t remaining_{0};
  std::uint64_t value_{0};
  bool has_value_{false};
};

/// Multiway-merge sorted `.kset` paths into one `.kuniv` (fan-in limited by group_size).
Error build_universe_from_presence_sets(const std::vector<std::string>& kset_paths,
                                        std::size_t kmer_size, const std::string& output_path,
                                        const std::string& tmpdir, std::size_t group_size = 32);

bool path_looks_universe(const std::string& path);

}  // namespace kmat
