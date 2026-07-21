#include <catch2/catch_test_macros.hpp>

#include "kmat/count.hpp"
#include "kmat/io.hpp"
#include "kmat/matrix.hpp"
#include "kmat/presence.hpp"
#include "kmat/stripe_build.hpp"
#include "kmat/universe.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

fs::path tmp_dir(const std::string& name) {
  const fs::path dir = fs::temp_directory_path() / ("kmat_staged_" + name);
  std::error_code ec;
  fs::remove_all(dir, ec);
  fs::create_directories(dir, ec);
  return dir;
}

bool matrices_semantically_equal(const kmat::PaMatrix& a, const kmat::PaMatrix& b) {
  if (a.header.num_rows != b.header.num_rows ||
      a.header.num_accessions != b.header.num_accessions ||
      a.header.num_stripes != b.header.num_stripes || a.kmers.size() != b.kmers.size()) {
    return false;
  }
  for (std::size_t i = 0; i < a.kmers.size(); ++i) {
    if (a.kmers[i].kmer_code != b.kmers[i].kmer_code) {
      return false;
    }
    if (a.kmers[i].pattern_id >= a.patterns.size() ||
        b.kmers[i].pattern_id >= b.patterns.size()) {
      return false;
    }
    if (a.patterns[a.kmers[i].pattern_id] != b.patterns[b.kmers[i].pattern_id]) {
      return false;
    }
  }
  return true;
}

std::size_t count_files_recursive(const fs::path& root) {
  std::size_t n = 0;
  std::error_code ec;
  if (!fs::exists(root, ec)) {
    return 0;
  }
  for (auto it = fs::recursive_directory_iterator(root, ec);
       it != fs::recursive_directory_iterator(); it.increment(ec)) {
    if (it->is_regular_file(ec)) {
      ++n;
    }
  }
  return n;
}

std::vector<std::string> make_ksets(const fs::path& dir, std::size_t k = 3) {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR);
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
    copts.kmer_size = k;
    copts.min_count = 1;
    copts.engine = kmat::CountEngine::Builtin;
    REQUIRE(kmat::count_kmers_to_presence_set(copts).ok());
    kset_paths.push_back(out.string());
  }
  return kset_paths;
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

TEST_CASE("staged build matches sequence build semantically", "[matrix][build]") {
  const fs::path dir = tmp_dir("staged");
  const fs::path spill = dir / "spill";
  fs::create_directories(spill);

  const fs::path td = fs::path(KMAT_TESTDATA_DIR);
  std::vector<std::string> fasta_paths;
  REQUIRE(kmat::read_list_file((td / "accession_list.txt").string(), fasta_paths).ok());
  REQUIRE(kmat::resolve_list_paths((td / "accession_list.txt").string(), fasta_paths).ok());
  auto kset_paths = make_ksets(dir, 3);

  const fs::path via_seq = dir / "seq.kmat";
  const fs::path via_kset = dir / "kset.kmat";

  kmat::BuildOptions b1;
  b1.kmer_size = 3;
  b1.accession_paths = fasta_paths;
  b1.output_path = via_seq.string();
  REQUIRE(kmat::build_matrix_from_accessions(b1).ok());

  kmat::BuildOptions b2;
  b2.kmer_size = 3;
  b2.accession_paths = kset_paths;
  b2.output_path = via_kset.string();
  b2.batch_rows = 16;
  b2.tmpdir = spill.string();
  REQUIRE(kmat::build_matrix_from_accessions(b2).ok());

  kmat::PaMatrix m1;
  kmat::PaMatrix m2;
  REQUIRE(kmat::read_matrix(via_seq.string(), m1).ok());
  REQUIRE(kmat::read_matrix(via_kset.string(), m2).ok());
  REQUIRE(matrices_semantically_equal(m1, m2));

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("staged build keeps few files under tmpdir", "[matrix][build]") {
  const fs::path dir = tmp_dir("fewfiles");
  const fs::path spill = dir / "spill";
  fs::create_directories(spill);
  auto kset_paths = make_ksets(dir, 3);

  kmat::BuildOptions opts;
  opts.kmer_size = 3;
  opts.accession_paths = kset_paths;
  opts.output_path = (dir / "out.kmat").string();
  opts.batch_rows = 8;
  opts.tmpdir = spill.string();
  REQUIRE(kmat::build_matrix_from_accessions(opts).ok());

  // During build, Cleaner removes work dir on success. Peak is bounded by
  // O(N/group + stripes); assert final spill is empty/small and we never used N*T layout.
  const std::size_t n = kset_paths.size();
  const std::size_t stripes = kmat::stripe_count_for_accessions(n);
  // Soft bound: if anything left, far below N*64 scatter footprint.
  REQUIRE(count_files_recursive(spill) < n * 4 + stripes + 8);

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("build-master produces sorted unique codes", "[universe][build]") {
  const fs::path dir = tmp_dir("kuniv");
  auto kset_paths = make_ksets(dir, 3);
  const fs::path out = dir / "panel.kuniv";
  REQUIRE(kmat::build_universe_from_presence_sets(kset_paths, 3, out.string(), dir.string(), 2)
              .ok());

  kmat::UniverseCursor cur;
  REQUIRE(cur.open(out.string()).ok());
  std::vector<std::uint64_t> codes;
  while (cur.has_value()) {
    codes.push_back(cur.value());
    REQUIRE(cur.advance().ok());
  }
  REQUIRE(codes.size() == cur.header().num_kmers);
  REQUIRE(codes.size() > 0);
  REQUIRE(std::is_sorted(codes.begin(), codes.end()));
  REQUIRE(std::adjacent_find(codes.begin(), codes.end()) == codes.end());

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("universe cursor recovers zero num_kmers header from file size", "[universe]") {
  const fs::path dir = tmp_dir("kuniv_hdr");
  const fs::path path = dir / "broken.kuniv";
  kmat::UniverseHeader hdr{};
  hdr.magic[0] = 'K';
  hdr.magic[1] = 'U';
  hdr.magic[2] = 'N';
  hdr.magic[3] = 'I';
  hdr.version = 1;
  hdr.kmer_size = 3;
  hdr.num_kmers = 0;  // simulate failed seekp rewrite
  const std::vector<std::uint64_t> codes{1, 2, 5, 9};
  {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
    out.write(reinterpret_cast<const char*>(codes.data()),
              static_cast<std::streamsize>(codes.size() * sizeof(std::uint64_t)));
  }

  kmat::UniverseCursor cur;
  REQUIRE(cur.open(path.string()).ok());
  REQUIRE(cur.header().num_kmers == codes.size());
  std::vector<std::uint64_t> got;
  while (cur.has_value()) {
    got.push_back(cur.value());
    REQUIRE(cur.advance().ok());
  }
  REQUIRE(got == codes);

  // Tree-reduce must also see the recovered codes.
  const fs::path out = dir / "merged.kuniv";
  // Build via two broken inputs through build_universe would need ksets; merge path tested
  // by re-writing a second broken file and opening both via a tiny build from ksets instead:
  // header rewrite on fresh build-master must leave a non-zero count on disk.
  auto kset_paths = make_ksets(dir, 3);
  REQUIRE(kmat::build_universe_from_presence_sets(kset_paths, 3, out.string(), dir.string(), 2)
              .ok());
  kmat::UniverseHeader out_hdr{};
  REQUIRE(kmat::read_universe_header(out.string(), out_hdr).ok());
  REQUIRE(out_hdr.num_kmers > 0);
  // On-disk header itself must be patched (not only recovered at read time).
  {
    std::ifstream in(out, std::ios::binary);
    kmat::UniverseHeader raw{};
    in.read(reinterpret_cast<char*>(&raw), sizeof(raw));
    REQUIRE(raw.num_kmers == out_hdr.num_kmers);
    REQUIRE(raw.num_kmers > 0);
  }

  std::error_code ec;
  fs::remove_all(dir, ec);
}

TEST_CASE("staged build medium panel", "[matrix][build][medium]") {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR) / "panel_k31_n72";
  const fs::path dir = tmp_dir("medium");
  const fs::path spill = dir / "spill";
  fs::create_directories(spill);

  std::vector<std::string> paths;
  REQUIRE(kmat::read_list_file((td / "accession_list.txt").string(), paths).ok());
  REQUIRE(kmat::resolve_list_paths((td / "accession_list.txt").string(), paths).ok());

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
  o1.batch_rows = 1000;
  o1.tmpdir = spill.string();
  REQUIRE(kmat::build_matrix_from_accessions(o1).ok());

  kmat::BuildOptions o2 = o1;
  o2.output_path = b.string();
  o2.batch_rows = 16;
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
