#include <catch2/catch_test_macros.hpp>

#include "kmat/io.hpp"
#include "kmat/kmer.hpp"
#include "kmat/matrix.hpp"

#include <cstdint>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

TEST_CASE("build_matrix_from_sequences writes v2 pattern matrix", "[matrix]") {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR);
  const fs::path list_path = td / "accession_list.txt";
  const fs::path out_path = fs::temp_directory_path() / "kmat_test_panel.kmat";

  std::vector<std::string> paths;
  REQUIRE(kmat::read_list_file(list_path.string(), paths).ok());
  REQUIRE(kmat::resolve_list_paths(list_path.string(), paths).ok());

  kmat::BuildOptions opts;
  opts.kmer_size = 3;
  opts.accession_paths = std::move(paths);
  opts.output_path = out_path.string();
  REQUIRE(kmat::build_matrix_from_sequences(opts).ok());

  kmat::PaMatrix matrix;
  REQUIRE(kmat::read_matrix(out_path.string(), matrix).ok());
  REQUIRE(matrix.header.version == 2);
  REQUIRE(matrix.header.kmer_size == 3);
  REQUIRE(matrix.header.num_accessions == 6);
  REQUIRE(matrix.header.num_stripes == 1);
  REQUIRE(matrix.kmers.size() > 0);
  REQUIRE(matrix.kmers.size() == matrix.header.num_rows);
  REQUIRE(matrix.patterns.size() == matrix.header.reserved);
  REQUIRE(matrix.patterns.size() > 0);
  REQUIRE(matrix.patterns.size() <= matrix.kmers.size());

  for (std::size_t i = 1; i < matrix.kmers.size(); ++i) {
    REQUIRE(matrix.kmers[i].kmer_code > matrix.kmers[i - 1].kmer_code);
  }

  std::error_code ec;
  fs::remove(out_path, ec);
}

TEST_CASE("presence_bitstring uses accession order", "[matrix]") {
  kmat::MatrixRow row;
  row.words = {0b000101ULL};  // bits 0 and 2 set
  REQUIRE(kmat::presence_bitstring(row, 3) == "101");
  REQUIRE(kmat::count_carriers(row, 3) == 2);
}

TEST_CASE("encode known k-mer round-trips through matrix row", "[matrix]") {
  std::uint64_t code = 0;
  REQUIRE(kmat::encode_kmer("ACG", code).ok());

  kmat::MatrixRow row;
  row.kmer_code = code;
  row.words = {7ULL};  // present in first three accessions
  REQUIRE(kmat::presence_bitstring(row, 6) == "111000");
}

TEST_CASE("medium panel has fewer patterns than k-mers", "[matrix][medium]") {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR) / "panel_k31_n72";
  const fs::path list_path = td / "accession_list.txt";
  const fs::path out_path = fs::temp_directory_path() / "kmat_test_panel_k31.kmat";

  std::vector<std::string> paths;
  REQUIRE(kmat::read_list_file(list_path.string(), paths).ok());
  REQUIRE(kmat::resolve_list_paths(list_path.string(), paths).ok());

  kmat::BuildOptions opts;
  opts.kmer_size = 31;
  opts.accession_paths = paths;
  opts.output_path = out_path.string();
  REQUIRE(kmat::build_matrix_from_sequences(opts).ok());

  kmat::PaMatrix matrix;
  REQUIRE(kmat::read_matrix(out_path.string(), matrix).ok());
  REQUIRE(matrix.header.version == 2);
  REQUIRE(matrix.header.num_stripes == 2);
  REQUIRE(matrix.patterns.size() < matrix.kmers.size());

  std::error_code ec;
  fs::remove(out_path, ec);
}
