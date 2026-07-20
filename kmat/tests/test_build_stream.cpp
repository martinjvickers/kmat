#include <catch2/catch_test_macros.hpp>

#include "kmat/count.hpp"
#include "kmat/io.hpp"
#include "kmat/matrix.hpp"
#include "kmat/presence.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

fs::path tmp_dir(const std::string& name) {
  const fs::path dir = fs::temp_directory_path() / ("kmat_build_stream_" + name);
  std::error_code ec;
  fs::remove_all(dir, ec);
  fs::create_directories(dir, ec);
  return dir;
}

bool matrices_semantically_equal(const kmat::PaMatrix& a, const kmat::PaMatrix& b) {
  if (a.header.num_rows != b.header.num_rows || a.header.reserved != b.header.reserved ||
      a.header.num_accessions != b.header.num_accessions ||
      a.header.num_stripes != b.header.num_stripes || a.kmers.size() != b.kmers.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.kmers.size(); ++i) {
    if (a.kmers[i].kmer_code != b.kmers[i].kmer_code) {
      return false;
    }
    if (a.patterns[a.kmers[i].pattern_id] != b.patterns[b.kmers[i].pattern_id]) {
      return false;
    }
  }
  return true;
}

}  // namespace

TEST_CASE("PresenceSetCursor streams sorted k-mers", "[presence][build]") {
  const fs::path dir = tmp_dir("cursor");
  const fs::path path = dir / "a.kset";

  kmat::PresenceSet set;
  set.header.kmer_size = 3;
  set.header.min_count = 1;
  set.kmers = {1, 5, 9, 20};
  set.header.num_kmers = set.kmers.size();
  REQUIRE(kmat::write_presence_set(path.string(), set).ok());

  kmat::PresenceSetCursor cur;
  REQUIRE(cur.open(path.string()).ok());
  std::vector<std::uint64_t> got;
  while (cur.has_value()) {
    got.push_back(cur.value());
    REQUIRE(cur.advance().ok());
  }
  REQUIRE(got == set.kmers);

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("streaming build matches under tiny memory budget", "[matrix][build]") {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR);
  const fs::path dir = tmp_dir("mem");

  std::vector<std::string> fasta_paths;
  REQUIRE(kmat::read_list_file((td / "accession_list.txt").string(), fasta_paths).ok());
  REQUIRE(kmat::resolve_list_paths((td / "accession_list.txt").string(), fasta_paths).ok());

  std::vector<std::string> kset_paths;
  for (const std::string& fa : fasta_paths) {
    const std::string id = kmat::accession_id_from_path(fa);
    const fs::path out = dir / (id + ".kset");
    kmat::CountOptions copts;
    copts.input_path = fa;
    copts.output_path = out.string();
    copts.kmer_size = 3;
    copts.min_count = 1;
    copts.engine = kmat::CountEngine::Builtin;
    REQUIRE(kmat::count_kmers_to_presence_set(copts).ok());
    kset_paths.push_back(out.string());
  }

  const fs::path large = dir / "large.kmat";
  const fs::path tiny = dir / "tiny.kmat";
  const fs::path spill = dir / "spill";
  fs::create_directories(spill);

  kmat::BuildOptions bloated;
  bloated.kmer_size = 3;
  bloated.accession_paths = kset_paths;
  bloated.output_path = large.string();
  bloated.memory_bytes = 8ull << 30;
  bloated.batch_rows = 100000;
  bloated.tmpdir = spill.string();
  bloated.num_threads = 2;
  REQUIRE(kmat::build_matrix_from_accessions(bloated).ok());

  kmat::BuildOptions tight;
  tight.kmer_size = 3;
  tight.accession_paths = kset_paths;
  tight.output_path = tiny.string();
  tight.memory_bytes = 1ull << 20;  // 1 MiB → many partitions
  tight.batch_rows = 2;
  tight.tmpdir = spill.string();
  tight.num_threads = 4;
  REQUIRE(kmat::build_matrix_from_accessions(tight).ok());

  kmat::PaMatrix m1;
  kmat::PaMatrix m2;
  REQUIRE(kmat::read_matrix(large.string(), m1).ok());
  REQUIRE(kmat::read_matrix(tiny.string(), m2).ok());
  REQUIRE(matrices_semantically_equal(m1, m2));

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("streaming build medium panel under constrained memory", "[matrix][build][medium]") {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR) / "panel_k31_n72";
  const fs::path dir = tmp_dir("medium_mem");
  const fs::path spill = dir / "spill";
  fs::create_directories(spill);

  std::vector<std::string> paths;
  REQUIRE(kmat::read_list_file((td / "accession_list.txt").string(), paths).ok());
  REQUIRE(kmat::resolve_list_paths((td / "accession_list.txt").string(), paths).ok());

  // Count a subset to .kset for a faster streaming check (first 8 accessions).
  const std::size_t n = std::min<std::size_t>(8, paths.size());
  std::vector<std::string> kset_paths;
  for (std::size_t i = 0; i < n; ++i) {
    const std::string id = kmat::accession_id_from_path(paths[i]);
    const fs::path out = dir / (id + ".kset");
    kmat::CountOptions copts;
    copts.input_path = paths[i];
    copts.output_path = out.string();
    copts.kmer_size = 31;
    copts.min_count = 1;
    copts.engine = kmat::CountEngine::Builtin;
    REQUIRE(kmat::count_kmers_to_presence_set(copts).ok());
    kset_paths.push_back(out.string());
  }

  const fs::path a = dir / "a.kmat";
  const fs::path b = dir / "b.kmat";

  kmat::BuildOptions o1;
  o1.kmer_size = 31;
  o1.accession_paths = kset_paths;
  o1.output_path = a.string();
  o1.memory_bytes = 512ull << 20;
  o1.batch_rows = 1000;
  o1.tmpdir = spill.string();
  REQUIRE(kmat::build_matrix_from_accessions(o1).ok());

  kmat::BuildOptions o2;
  o2.kmer_size = 31;
  o2.accession_paths = kset_paths;
  o2.output_path = b.string();
  o2.memory_bytes = 4ull << 20;
  o2.batch_rows = 16;
  o2.tmpdir = spill.string();
  o2.num_threads = 2;
  REQUIRE(kmat::build_matrix_from_accessions(o2).ok());

  kmat::PaMatrix m1;
  kmat::PaMatrix m2;
  REQUIRE(kmat::read_matrix(a.string(), m1).ok());
  REQUIRE(kmat::read_matrix(b.string(), m2).ok());
  REQUIRE(matrices_semantically_equal(m1, m2));
  REQUIRE(m1.header.num_accessions == n);

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("hierarchical merge matches direct merge under tiny FD budget", "[matrix][build]") {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR);
  const fs::path dir = tmp_dir("hier");
  const fs::path spill = dir / "spill";
  fs::create_directories(spill);

  std::vector<std::string> fasta_paths;
  REQUIRE(kmat::read_list_file((td / "accession_list.txt").string(), fasta_paths).ok());
  REQUIRE(kmat::resolve_list_paths((td / "accession_list.txt").string(), fasta_paths).ok());

  std::vector<std::string> kset_paths;
  for (const std::string& fa : fasta_paths) {
    const std::string id = kmat::accession_id_from_path(fa);
    const fs::path out = dir / (id + ".kset");
    kmat::CountOptions copts;
    copts.input_path = fa;
    copts.output_path = out.string();
    copts.kmer_size = 3;
    copts.min_count = 1;
    copts.engine = kmat::CountEngine::Builtin;
    REQUIRE(kmat::count_kmers_to_presence_set(copts).ok());
    kset_paths.push_back(out.string());
  }

  const fs::path direct = dir / "direct.kmat";
  const fs::path hier = dir / "hier.kmat";

  kmat::BuildOptions o1;
  o1.kmer_size = 3;
  o1.accession_paths = kset_paths;
  o1.output_path = direct.string();
  o1.memory_bytes = 64ull << 20;
  o1.batch_rows = 16;
  o1.tmpdir = spill.string();
  o1.num_threads = 2;
  unsetenv("KMAT_BUILD_MAX_OPEN");
  REQUIRE(kmat::build_matrix_from_accessions(o1).ok());

  // Force hierarchical accession merge (group size 2) regardless of ulimit.
  REQUIRE(setenv("KMAT_BUILD_MAX_OPEN", "2", 1) == 0);
  kmat::BuildOptions o2 = o1;
  o2.output_path = hier.string();
  REQUIRE(kmat::build_matrix_from_accessions(o2).ok());
  unsetenv("KMAT_BUILD_MAX_OPEN");

  kmat::PaMatrix m1;
  kmat::PaMatrix m2;
  REQUIRE(kmat::read_matrix(direct.string(), m1).ok());
  REQUIRE(kmat::read_matrix(hier.string(), m2).ok());
  REQUIRE(matrices_semantically_equal(m1, m2));

  std::error_code ec;
  fs::remove_all(dir, ec);
}
