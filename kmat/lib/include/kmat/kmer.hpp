#pragma once

#include "kmat/error.hpp"

#include <cstdint>
#include <functional>
#include <string>
#include <string_view>

namespace kmat {

/// LSB-first 2-bit encoding: A=00, C=01, G=10, T=11 (position i uses bits 2*i .. 2*i+1).
Error encode_kmer(std::string_view kmer, std::uint64_t& out);

Error decode_kmer(std::uint64_t code, std::size_t kmer_length, std::string& out);

/// Sliding-window encode with O(1) roll when possible. Skips windows containing invalid bases.
/// Invokes `fn(code)` for each valid k-mer (forward strand only).
Error for_each_encoded_kmer(std::string_view sequence, std::size_t kmer_size,
                            const std::function<void(std::uint64_t)>& fn);

}  // namespace kmat
