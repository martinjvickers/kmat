#pragma once

#include "kmat/error.hpp"

#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace kmat {

Error reverse_complement(std::string_view sequence, std::string& out);

/// Extract forward-strand k-mers from a DNA sequence (min count 1).
Error kmer_set_from_sequence(std::string_view sequence, std::size_t kmer_size,
                             std::unordered_set<std::uint64_t>& out);

/// Sliding-window k-mers from a sequence plus reverse complements.
Error gene_kmer_set(std::string_view sequence, std::size_t kmer_size,
                    std::unordered_set<std::uint64_t>& out);

}  // namespace kmat
