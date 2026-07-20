#include <catch2/catch_test_macros.hpp>

#include "kmat/io.hpp"
#include "kmat/kmer.hpp"

#include <cstdint>
#include <filesystem>
#include <string>

namespace {

constexpr std::uint64_t kExpectedAcgt = 228ULL;  // A=00,C=01,G=10,T=11 at positions 0..3

}  // namespace

TEST_CASE("encode_kmer packs LSB-first 2-bit bases", "[kmer]") {
  std::uint64_t code = 0;
  const auto err = kmat::encode_kmer("ACGT", code);
  REQUIRE(err.ok());
  REQUIRE(code == kExpectedAcgt);
}

TEST_CASE("decode_kmer inverts encode_kmer", "[kmer]") {
  std::uint64_t code = 0;
  REQUIRE(kmat::encode_kmer("ACGT", code).ok());

  std::string decoded;
  REQUIRE(kmat::decode_kmer(code, 4, decoded).ok());
  REQUIRE(decoded == "ACGT");
}

TEST_CASE("encode_kmer rejects invalid nucleotides", "[kmer]") {
  std::uint64_t code = 0;
  const auto err = kmat::encode_kmer("ACGX", code);
  REQUIRE(!err.ok());
  REQUIRE(err.code == kmat::ErrorCode::InvalidKmer);
}

TEST_CASE("encode/decode round-trip for 31-mer", "[kmer]") {
  std::string kmer;
  kmer.reserve(31);
  for (std::size_t i = 0; i < 31; ++i) {
    switch (i % 4) {
      case 0:
        kmer.push_back('A');
        break;
      case 1:
        kmer.push_back('C');
        break;
      case 2:
        kmer.push_back('G');
        break;
      case 3:
        kmer.push_back('T');
        break;
    }
  }

  std::uint64_t code = 0;
  REQUIRE(kmat::encode_kmer(kmer, code).ok());

  std::string out;
  REQUIRE(kmat::decode_kmer(code, kmer.size(), out).ok());
  REQUIRE(out == kmer);
}

TEST_CASE("read_list_file skips blanks and comments", "[io]") {
  const std::filesystem::path fixture =
      std::filesystem::path(KMAT_TESTDATA_DIR) / "accession_list.txt";

  std::vector<std::string> lines;
  const auto err = kmat::read_list_file(fixture.string(), lines);
  REQUIRE(err.ok());
  REQUIRE(lines.size() == 6);
  REQUIRE(lines[0] == "accessions/sample_a.fasta");
  REQUIRE(lines[5] == "accessions/sample_f.fasta");
}
