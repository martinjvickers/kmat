#include "kmat/kmer.hpp"

namespace kmat {

namespace {

constexpr int base_bits(char nucleotide, bool& ok) {
  switch (nucleotide) {
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
    default:
      ok = false;
      return 0;
  }
}

constexpr char bits_to_base(unsigned bits) {
  switch (bits & 3U) {
    case 0:
      return 'A';
    case 1:
      return 'C';
    case 2:
      return 'G';
    case 3:
      return 'T';
    default:
      return '?';
  }
}

}  // namespace

Error encode_kmer(std::string_view kmer, std::uint64_t& out) {
  if (kmer.empty()) {
    return Error::invalid_argument("k-mer length must be greater than zero");
  }
  if (kmer.size() > 32) {
    return Error::invalid_argument("k-mer length exceeds 64-bit packing limit (32)");
  }

  out = 0;
  for (std::size_t i = 0; i < kmer.size(); ++i) {
    bool ok = true;
    const int bits = base_bits(kmer[i], ok);
    if (!ok) {
      return Error::invalid_kmer(std::string("invalid nucleotide in k-mer: ") + kmer[i]);
    }
    out |= static_cast<std::uint64_t>(bits) << (2 * i);
  }
  return Error::success();
}

Error decode_kmer(std::uint64_t code, std::size_t kmer_length, std::string& out) {
  if (kmer_length == 0) {
    return Error::invalid_argument("k-mer length must be greater than zero");
  }
  if (kmer_length > 32) {
    return Error::invalid_argument("k-mer length exceeds 64-bit packing limit (32)");
  }

  out.clear();
  out.reserve(kmer_length);
  for (std::size_t i = 0; i < kmer_length; ++i) {
    const unsigned bits = static_cast<unsigned>((code >> (2 * i)) & 3ULL);
    out.push_back(bits_to_base(bits));
  }
  return Error::success();
}

Error for_each_encoded_kmer(std::string_view sequence, std::size_t kmer_size,
                            const std::function<void(std::uint64_t)>& fn) {
  if (kmer_size == 0 || kmer_size > 32) {
    return Error::invalid_argument("k-mer size must be in 1..32");
  }
  if (sequence.size() < kmer_size) {
    return Error::success();
  }

  const std::uint64_t high_shift = 2ULL * (kmer_size - 1);
  std::uint64_t code = 0;
  std::size_t run = 0;

  for (std::size_t i = 0; i < sequence.size(); ++i) {
    bool ok = true;
    const int bits = base_bits(sequence[i], ok);
    if (!ok) {
      run = 0;
      code = 0;
      continue;
    }

    if (run < kmer_size) {
      code |= static_cast<std::uint64_t>(bits) << (2 * run);
      ++run;
      if (run == kmer_size) {
        fn(code);
      }
    } else {
      // Drop oldest base (LSB pair); append new base at the high end.
      code >>= 2;
      code |= static_cast<std::uint64_t>(bits) << high_shift;
      fn(code);
    }
  }
  return Error::success();
}

}  // namespace kmat
