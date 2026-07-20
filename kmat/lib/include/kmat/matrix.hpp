#pragma once

#include "kmat/error.hpp"
#include "kmat/matrix_layout.hpp"

#include <cstdint>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

namespace kmat {

/// On-disk header (40 bytes), little-endian. See docs/FORMAT.md.
struct MatrixHeader {
  char magic[4]{'K', 'M', 'A', 'T'};
  std::uint32_t version{2};
  std::uint32_t kmer_size{0};
  std::uint32_t num_accessions{0};
  std::uint32_t num_stripes{0};
  /// v1: k-mer row count. v2: k-mer map entry count.
  std::uint64_t num_rows{0};
  /// v1: unused. v2: number of unique PA patterns.
  std::uint64_t reserved{0};
};

static_assert(sizeof(MatrixHeader) == 40, "MatrixHeader must be 40 bytes");

#pragma pack(push, 1)
struct KmerMapEntry {
  std::uint64_t kmer_code{0};
  std::uint32_t pattern_id{0};
  std::uint32_t pad{0};
};
#pragma pack(pop)

static_assert(sizeof(KmerMapEntry) == 16, "KmerMapEntry must be 16 bytes");

/// Dense row view (v1 legacy / helpers).
struct MatrixRow {
  std::uint64_t kmer_code{0};
  std::vector<std::uint64_t> words;
};

/// Unified in-memory matrix (always pattern + k-mer map after load).
struct PaMatrix {
  MatrixHeader header{};
  std::vector<std::vector<std::uint64_t>> patterns;
  std::vector<KmerMapEntry> kmers;

  std::size_t num_patterns() const { return patterns.size(); }
  std::size_t num_kmers() const { return kmers.size(); }
};

struct BuildOptions {
  std::size_t kmer_size{31};
  std::vector<std::string> accession_paths;
  std::string output_path;
  /// 0 = use process runtime profile thread count.
  std::size_t num_threads{0};
};

Error build_matrix_from_sequences(const BuildOptions& opts);

/// Build v2 matrix from per-accession `.kset` presence sets (same list order as GWAS).
Error build_matrix_from_presence_sets(const BuildOptions& opts);

/// Dispatch: all `.kset` → presence build; otherwise sequence files (FASTA/FASTQ/.gz).
Error build_matrix_from_accessions(const BuildOptions& opts);

/// Alias for build_matrix_from_sequences.
Error build_matrix_from_fastas(const BuildOptions& opts);

Error read_matrix(const std::string& path, PaMatrix& matrix);

Error write_matrix(const std::string& path, const PaMatrix& matrix);

Error read_matrix_list(const std::string& list_path, std::vector<std::string>& stripe_paths);

/// Concatenate v1 stripe files into one in-memory matrix (inflated to patterns).
Error load_matrix_from_list(const std::string& matrix_list_path, std::size_t num_accessions,
                            PaMatrix& matrix);

bool get_presence_bit(const std::vector<std::uint64_t>& words, std::size_t accession_index);

bool get_presence_bit(const MatrixRow& row, std::size_t accession_index);

std::string presence_bitstring(const std::vector<std::uint64_t>& words, std::size_t num_accessions);

std::string presence_bitstring(const MatrixRow& row, std::size_t num_accessions);

std::size_t count_carriers(const std::vector<std::uint64_t>& words, std::size_t num_accessions);

std::size_t count_carriers(const MatrixRow& row, std::size_t num_accessions);

/// Binary search k-mer map; returns pattern_id if found.
std::optional<std::uint32_t> find_pattern_id(const PaMatrix& matrix, std::uint64_t kmer_code);

/// Collect all k-mer codes that share a pattern_id (matrix.kmers must be sorted by code; scan is O(K)).
std::vector<std::uint64_t> kmers_for_pattern(const PaMatrix& matrix, std::uint32_t pattern_id);

/// Build inverted index pattern_id → k-mer codes (O(K) once).
std::vector<std::vector<std::uint64_t>> pattern_kmer_index(const PaMatrix& matrix);

}  // namespace kmat
