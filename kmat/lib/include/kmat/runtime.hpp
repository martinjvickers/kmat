#pragma once

#include <cstddef>
#include <string>

namespace kmat {

enum class RuntimeProfile {
  Laptop,
  Hpc,
};

struct RuntimeConfig {
  RuntimeProfile profile{RuntimeProfile::Laptop};
  /// Worker threads for CPU-bound loops (GWAS patterns, build ingest). 0 = resolve from profile.
  std::size_t num_threads{0};
  /// Preferred iostream buffer size for large sequential matrix I/O.
  std::size_t io_buffer_bytes{1u << 20};
};

/// Resolve profile defaults (threads / buffer). Explicit num_threads > 0 wins.
RuntimeConfig resolve_runtime(RuntimeProfile profile, std::size_t num_threads_override = 0);

/// Process-wide defaults (CLI `--profile` / `--threads`).
void set_runtime_config(const RuntimeConfig& cfg);
RuntimeConfig runtime_config();

RuntimeProfile parse_runtime_profile(const std::string& name, bool& ok);
const char* runtime_profile_name(RuntimeProfile profile);

std::size_t effective_threads(const RuntimeConfig& cfg);

/// Split [begin, end) across threads; fn(i) called for each index. Deterministic work partition.
template <typename Fn>
void parallel_for(std::size_t begin, std::size_t end, std::size_t num_threads, Fn&& fn);

}  // namespace kmat

#include "kmat/runtime_inl.hpp"
