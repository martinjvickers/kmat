#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "kmat/gwas.hpp"
#include "kmat/io.hpp"
#include "kmat/kmer.hpp"
#include "kmat/matrix.hpp"
#include "kmat/pca.hpp"
#include "kmat/runtime.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

TEST_CASE("resolve_runtime laptop vs hpc", "[runtime]") {
  const auto laptop = kmat::resolve_runtime(kmat::RuntimeProfile::Laptop, 0);
  const auto hpc = kmat::resolve_runtime(kmat::RuntimeProfile::Hpc, 0);
  REQUIRE(laptop.io_buffer_bytes < hpc.io_buffer_bytes);
  REQUIRE(laptop.num_threads >= 1);
  REQUIRE(hpc.num_threads >= laptop.num_threads);

  const auto forced = kmat::resolve_runtime(kmat::RuntimeProfile::Hpc, 3);
  REQUIRE(forced.num_threads == 3);
}

TEST_CASE("parse_runtime_profile accepts aliases", "[runtime]") {
  bool ok = false;
  REQUIRE(kmat::parse_runtime_profile("laptop", ok) == kmat::RuntimeProfile::Laptop);
  REQUIRE(ok);
  REQUIRE(kmat::parse_runtime_profile("HPC", ok) == kmat::RuntimeProfile::Hpc);
  REQUIRE(ok);
  REQUIRE(kmat::parse_runtime_profile("nope", ok) == kmat::RuntimeProfile::Laptop);
  REQUIRE_FALSE(ok);
}

TEST_CASE("for_each_encoded_kmer matches encode_kmer windows", "[kmer]") {
  const std::string seq = "ACGTNACGT";
  std::vector<std::uint64_t> rolled;
  REQUIRE(kmat::for_each_encoded_kmer(seq, 3, [&](std::uint64_t code) { rolled.push_back(code); })
              .ok());

  std::vector<std::uint64_t> naive;
  for (std::size_t i = 0; i + 3 <= seq.size(); ++i) {
    std::uint64_t code = 0;
    if (kmat::encode_kmer(seq.substr(i, 3), code).ok()) {
      naive.push_back(code);
    }
  }
  REQUIRE(rolled == naive);
}

TEST_CASE("parallel GWAS matches single-thread results", "[gwas][runtime]") {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR);
  const fs::path dir = fs::temp_directory_path() / "kmat_gwas_par";
  fs::remove_all(dir);
  fs::create_directories(dir);
  const fs::path matrix_path = dir / "panel.kmat";
  const fs::path pop_path = dir / "pop.tsv";

  std::vector<std::string> paths;
  REQUIRE(kmat::read_list_file((td / "accession_list.txt").string(), paths).ok());
  REQUIRE(kmat::resolve_list_paths((td / "accession_list.txt").string(), paths).ok());

  kmat::BuildOptions bopts;
  bopts.kmer_size = 3;
  bopts.accession_paths = paths;
  bopts.output_path = matrix_path.string();
  bopts.num_threads = 1;
  REQUIRE(kmat::build_matrix_from_accessions(bopts).ok());

  kmat::PcaOptions popts;
  popts.matrix_path = matrix_path.string();
  popts.accession_list_path = (td / "accession_list.txt").string();
  popts.num_pcs = 0;  // all PCs
  kmat::PcaResult pop;
  REQUIRE(kmat::run_pca(popts, pop).ok());
  REQUIRE(pop.pcs.front().size() == paths.size());
  REQUIRE(kmat::write_pca_tsv(pop_path.string(), pop).ok());

  auto run = [&](std::size_t threads) {
    kmat::GwasOptions opts;
    opts.matrix_path = matrix_path.string();
    opts.accession_list_path = (td / "accession_list.txt").string();
    opts.phenotype_path = (td / "phenotypes.tsv").string();
    opts.pop_path = pop_path.string();
    opts.kmer_size = 3;
    opts.print_all = true;
    opts.num_threads = threads;
    kmat::GwasResult result;
    REQUIRE(kmat::run_gwas(opts, result).ok());
    return result;
  };

  const kmat::GwasResult one = run(1);
  const kmat::GwasResult many = run(4);
  REQUIRE(one.hits.size() == many.hits.size());
  for (std::size_t i = 0; i < one.hits.size(); ++i) {
    REQUIRE(one.hits[i].kmer_code == many.hits[i].kmer_code);
    REQUIRE(one.hits[i].num_carriers == many.hits[i].num_carriers);
    REQUIRE(one.hits[i].p_values.size() == many.hits[i].p_values.size());
    for (std::size_t j = 0; j < one.hits[i].p_values.size(); ++j) {
      REQUIRE(one.hits[i].p_values[j] == Catch::Approx(many.hits[i].p_values[j]).margin(1e-12));
    }
  }
}
