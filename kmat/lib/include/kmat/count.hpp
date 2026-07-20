#pragma once

#include "kmat/error.hpp"
#include "kmat/presence.hpp"

#include <cstddef>
#include <cstdint>
#include <string>

namespace kmat {

enum class CountEngine {
  /// Disk-backed KMC CLI (production default).
  Kmc,
  /// In-process hashmap (small fixtures / no KMC installed).
  Builtin,
};

struct CountOptions {
  std::string input_path;
  std::string output_path;
  std::size_t kmer_size{31};
  std::uint32_t min_count{1};  // KMC-style -ci
  CountEngine engine{CountEngine::Kmc};
  /// 0 → use process runtime profile thread count (KMC -t).
  std::size_t num_threads{0};
  /// Scratch for KMC DB (empty → $TMPDIR or /tmp).
  std::string tmpdir;
  /// Override binary names/paths (empty → search PATH for "kmc" / "kmc_tools").
  std::string kmc_bin;
  std::string kmc_tools_bin;
};

/// Count k-mers from FASTQ/FASTA (.gz ok); keep codes with count >= min_count; write `.kset`.
Error count_kmers_to_presence_set(const CountOptions& opts);

struct ImportKmersOptions {
  std::string input_path;  // one DNA k-mer per line (or hex uint64 with 0x prefix)
  std::string output_path;
  std::size_t kmer_size{31};
};

/// Migration helper without linking KMC: import a text k-mer list into a `.kset`.
Error import_kmers_text_to_presence_set(const ImportKmersOptions& opts);

/// True if `name` is executable on PATH (or is an absolute/relative path to an executable).
bool find_executable(const std::string& name, std::string& resolved_path);

}  // namespace kmat
