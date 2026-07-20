#include <catch2/catch_test_macros.hpp>

#include "kmat/fastq.hpp"
#include "kmat/io.hpp"
#include "kmat/sequence.hpp"

#include <filesystem>
#include <string>
#include <unordered_set>

namespace fs = std::filesystem;

TEST_CASE("accession_id_from_path strips compound extensions", "[io]") {
  REQUIRE(kmat::accession_id_from_path("accessions/acc_000.fastq.gz") == "acc_000");
  REQUIRE(kmat::accession_id_from_path("/tmp/foo.fa.gz") == "foo");
  REQUIRE(kmat::accession_id_from_path("bar.fasta") == "bar");
  REQUIRE(kmat::accession_id_from_path("acc_001.kset") == "acc_001");
}

TEST_CASE("read gzipped FASTQ from medium panel", "[fastq]") {
  const fs::path path =
      fs::path(KMAT_TESTDATA_DIR) / "panel_k31_n72" / "accessions" / "acc_000.fastq.gz";
  REQUIRE(fs::exists(path));

  std::vector<kmat::FastqRecord> records;
  REQUIRE(kmat::read_fastq_sequences(path.string(), records).ok());
  REQUIRE(records.size() == 1);
  REQUIRE(records.front().sequence.size() > 31);

  std::unordered_set<std::uint64_t> kmers;
  REQUIRE(kmat::kmer_set_from_sequence(records.front().sequence, 31, kmers).ok());
  REQUIRE(kmers.size() > 0);
}

TEST_CASE("read_sequence_file dispatches FASTQ.gz", "[fastq]") {
  const fs::path path =
      fs::path(KMAT_TESTDATA_DIR) / "panel_k31_n72" / "accessions" / "acc_001.fastq.gz";
  std::vector<kmat::FastaRecord> records;
  REQUIRE(kmat::read_sequence_file(path.string(), records).ok());
  REQUIRE(records.size() == 1);
  REQUIRE(!records.front().sequence.empty());
}
