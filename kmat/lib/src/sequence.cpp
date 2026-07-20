#include "kmat/sequence.hpp"

#include "kmat/kmer.hpp"

namespace kmat {

Error reverse_complement(std::string_view sequence, std::string& out) {
  out.clear();
  out.reserve(sequence.size());
  for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
    switch (*it) {
      case 'A':
      case 'a':
        out.push_back('T');
        break;
      case 'C':
      case 'c':
        out.push_back('G');
        break;
      case 'G':
      case 'g':
        out.push_back('C');
        break;
      case 'T':
      case 't':
        out.push_back('A');
        break;
      case 'N':
      case 'n':
        out.push_back('N');
        break;
      default:
        return Error::invalid_kmer(std::string("invalid nucleotide: ") + *it);
    }
  }
  return Error::success();
}

Error kmer_set_from_sequence(std::string_view sequence, std::size_t kmer_size,
                             std::unordered_set<std::uint64_t>& out) {
  return for_each_encoded_kmer(sequence, kmer_size, [&](std::uint64_t code) { out.insert(code); });
}

Error gene_kmer_set(std::string_view sequence, std::size_t kmer_size,
                    std::unordered_set<std::uint64_t>& out) {
  if (auto err = kmer_set_from_sequence(sequence, kmer_size, out); !err.ok()) {
    return err;
  }
  std::string rc;
  if (auto err = reverse_complement(sequence, rc); !err.ok()) {
    return err;
  }
  return kmer_set_from_sequence(rc, kmer_size, out);
}

}  // namespace kmat
