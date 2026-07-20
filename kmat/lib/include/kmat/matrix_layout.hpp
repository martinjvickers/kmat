#pragma once

#include <cstddef>
#include <cstdint>

namespace kmat {

/// On-disk stripe record layout (legacy-compatible production shape).
///
/// Each stripe file holds one presence word per master k-mer row when the stripe
/// lists at most kAccessionsPerStripe accessions. Rows are fixed-size records:
///
///   { uint64_t kmer_code; uint64_t presence_word; }  // 16 bytes total
///
/// Global accession index g maps to stripe g / 64 and local bit g % 64 (LSB-first
/// within the stripe word). Downstream tools concatenate stripe words in list order
/// to reconstruct the full presence vector across unified_list.txt.
inline constexpr std::size_t kAccessionsPerStripe = 64;
inline constexpr std::size_t kStripeRecordBytes = sizeof(std::uint64_t) * 2;
inline constexpr std::size_t kKmerCodeFieldBytes = sizeof(std::uint64_t);
inline constexpr std::size_t kPresenceWordFieldBytes = sizeof(std::uint64_t);

inline constexpr std::size_t stripe_count_for_accessions(std::size_t accession_count) {
  return (accession_count + kAccessionsPerStripe - 1) / kAccessionsPerStripe;
}

inline constexpr std::size_t words_per_row(std::size_t accession_count) {
  return stripe_count_for_accessions(accession_count);
}

}  // namespace kmat
