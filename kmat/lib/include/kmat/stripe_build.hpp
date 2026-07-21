#pragma once

#include "kmat/error.hpp"
#include "kmat/matrix.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace kmat {

inline constexpr std::size_t kDefaultBatchRows = 100000;

struct CreateStripeOptions {
  std::string universe_path;
  std::string output_path;
  std::size_t kmer_size{31};
  /// Global panel accession count (header field).
  std::size_t num_accessions{0};
  /// Stripe index in [0, ceil(N/64)).
  std::size_t stripe_index{0};
  std::size_t batch_rows{kDefaultBatchRows};
};

struct FillStripeOptions {
  std::string stripe_path;
  std::string kset_path;
  /// Local bit index within the stripe word (0..63).
  std::size_t local_index{0};
  std::size_t batch_rows{kDefaultBatchRows};
};

struct CompressOptions {
  std::vector<std::string> stripe_paths;
  std::size_t num_accessions{0};
  std::size_t kmer_size{31};
  std::string output_path;
  std::size_t memory_bytes{0};
  std::size_t batch_rows{kDefaultBatchRows};
  std::string tmpdir;
};

/// Write a blank v1 single-word stripe from `.kuniv` (legacy-compatible layout).
Error create_blank_stripe(const CreateStripeOptions& opts);

/// Set bit `local_index` for each master k-mer present in the accession `.kset`.
Error fill_stripe_column(const FillStripeOptions& opts);

/// Stream filled stripes → globally pattern-deduped v2 `.kmat`.
Error compress_stripes_to_v2(const CompressOptions& opts);

/// Single-node driver: master → create all stripes → fill all columns → compress.
/// Uses few files under tmpdir (O(N/64) stripes + 1 kuniv), not N×T spills.
Error build_matrix_from_presence_sets_staged(const BuildOptions& opts);

}  // namespace kmat
