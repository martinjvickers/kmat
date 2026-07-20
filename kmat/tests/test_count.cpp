#include <catch2/catch_test_macros.hpp>

#include "kmat/count.hpp"
#include "kmat/io.hpp"
#include "kmat/matrix.hpp"
#include "kmat/presence.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

fs::path tmp_dir(const std::string& name) {
  const fs::path dir = fs::temp_directory_path() / ("kmat_" + name);
  fs::remove_all(dir);
  fs::create_directories(dir);
  return dir;
}

}  // namespace

TEST_CASE("presence set round-trips", "[presence]") {
  const fs::path dir = tmp_dir("presence_rt");
  const fs::path path = dir / "acc.kset";

  kmat::PresenceSet set;
  set.header.kmer_size = 3;
  set.header.min_count = 2;
  set.kmers = {1ULL, 5ULL, 9ULL};
  set.header.num_kmers = set.kmers.size();
  REQUIRE(kmat::write_presence_set(path.string(), set).ok());

  kmat::PresenceSet loaded;
  REQUIRE(kmat::read_presence_set(path.string(), loaded).ok());
  REQUIRE(loaded.header.kmer_size == 3);
  REQUIRE(loaded.header.min_count == 2);
  REQUIRE(loaded.kmers == set.kmers);
}

TEST_CASE("count --ci filters low-abundance k-mers", "[count]") {
  const fs::path dir = tmp_dir("count_ci");
  const fs::path fasta = dir / "acc.fa";
  {
    std::ofstream out(fasta);
    // AAA appears twice (windows 0 and overlapping? AAA only length 3 once per AAA)
    // Use: AAACAAA → AAA at 0 and 4
    out << ">acc\nAAACAAA\n";
  }

  const fs::path kset1 = dir / "ci1.kset";
  kmat::CountOptions opts;
  opts.input_path = fasta.string();
  opts.output_path = kset1.string();
  opts.kmer_size = 3;
  opts.min_count = 1;
  REQUIRE(kmat::count_kmers_to_presence_set(opts).ok());

  kmat::PresenceSet all;
  REQUIRE(kmat::read_presence_set(kset1.string(), all).ok());
  REQUIRE(all.kmers.size() >= 2);

  const fs::path kset2 = dir / "ci2.kset";
  opts.output_path = kset2.string();
  opts.min_count = 2;
  REQUIRE(kmat::count_kmers_to_presence_set(opts).ok());

  kmat::PresenceSet filtered;
  REQUIRE(kmat::read_presence_set(kset2.string(), filtered).ok());
  REQUIRE(filtered.kmers.size() < all.kmers.size());
  REQUIRE(filtered.header.min_count == 2);
}

TEST_CASE("import-kmers text to presence set", "[count]") {
  const fs::path dir = tmp_dir("import_kmers");
  const fs::path txt = dir / "kmers.txt";
  {
    std::ofstream out(txt);
    out << "AAA\n";
    out << "AAC\n";
    out << "# comment\n";
    out << "AAA\n";  // duplicate
  }

  const fs::path kset = dir / "out.kset";
  kmat::ImportKmersOptions opts;
  opts.input_path = txt.string();
  opts.output_path = kset.string();
  opts.kmer_size = 3;
  REQUIRE(kmat::import_kmers_text_to_presence_set(opts).ok());

  kmat::PresenceSet set;
  REQUIRE(kmat::read_presence_set(kset.string(), set).ok());
  REQUIRE(set.kmers.size() == 2);
}

TEST_CASE("count then build matches direct sequence build at --ci 1", "[count][matrix]") {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR);
  const fs::path dir = tmp_dir("count_build");

  std::vector<std::string> fasta_paths;
  REQUIRE(kmat::read_list_file((td / "accession_list.txt").string(), fasta_paths).ok());
  REQUIRE(kmat::resolve_list_paths((td / "accession_list.txt").string(), fasta_paths).ok());

  std::vector<std::string> kset_paths;
  kset_paths.reserve(fasta_paths.size());
  for (const std::string& fa : fasta_paths) {
    const std::string id = kmat::accession_id_from_path(fa);
    const fs::path out = dir / (id + ".kset");
    kmat::CountOptions copts;
    copts.input_path = fa;
    copts.output_path = out.string();
    copts.kmer_size = 3;
    copts.min_count = 1;
    REQUIRE(kmat::count_kmers_to_presence_set(copts).ok());
    kset_paths.push_back(out.string());
  }

  const fs::path direct = dir / "direct.kmat";
  const fs::path via_count = dir / "via_count.kmat";

  kmat::BuildOptions b1;
  b1.kmer_size = 3;
  b1.accession_paths = fasta_paths;
  b1.output_path = direct.string();
  REQUIRE(kmat::build_matrix_from_accessions(b1).ok());

  kmat::BuildOptions b2;
  b2.kmer_size = 3;
  b2.accession_paths = kset_paths;
  b2.output_path = via_count.string();
  REQUIRE(kmat::build_matrix_from_accessions(b2).ok());

  kmat::PaMatrix m1;
  kmat::PaMatrix m2;
  REQUIRE(kmat::read_matrix(direct.string(), m1).ok());
  REQUIRE(kmat::read_matrix(via_count.string(), m2).ok());
  REQUIRE(m1.header.num_rows == m2.header.num_rows);
  REQUIRE(m1.header.reserved == m2.header.reserved);
  REQUIRE(m1.kmers.size() == m2.kmers.size());
  for (std::size_t i = 0; i < m1.kmers.size(); ++i) {
    REQUIRE(m1.kmers[i].kmer_code == m2.kmers[i].kmer_code);
    REQUIRE(m1.patterns[m1.kmers[i].pattern_id] == m2.patterns[m2.kmers[i].pattern_id]);
  }
}
