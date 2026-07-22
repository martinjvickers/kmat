#include <catch2/catch_test_macros.hpp>

#include "kmat/gene.hpp"
#include "kmat/gwas.hpp"
#include "kmat/io.hpp"
#include "kmat/matrix.hpp"
#include "kmat/pca.hpp"
#include "kmat/validate.hpp"

#include <filesystem>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

namespace {

fs::path temp_dir_medium() {
  const fs::path base = fs::temp_directory_path() / "kmat_integration_k31_n72";
  std::error_code ec;
  fs::create_directories(base, ec);
  return base;
}

}  // namespace

TEST_CASE("end-to-end build pop gwas gene on k31 n72 FASTQ.gz panel", "[integration][medium]") {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR) / "panel_k31_n72";
  const fs::path out = temp_dir_medium();

  std::vector<std::string> paths;
  const std::string list_path = (td / "accession_list.txt").string();
  REQUIRE(kmat::read_list_file(list_path, paths).ok());
  REQUIRE(kmat::resolve_list_paths(list_path, paths).ok());
  REQUIRE(paths.size() == 72);

  const std::string matrix_path = (out / "panel.kmat").string();
  kmat::BuildOptions build_opts;
  build_opts.kmer_size = 31;
  build_opts.accession_paths = paths;
  build_opts.output_path = matrix_path;
  REQUIRE(kmat::build_matrix_from_sequences(build_opts).ok());

  kmat::PaMatrix matrix;
  REQUIRE(kmat::read_matrix(matrix_path, matrix).ok());
  REQUIRE(matrix.header.num_accessions == 72);
  REQUIRE(matrix.header.num_stripes == 2);
  REQUIRE(matrix.header.kmer_size == 31);
  REQUIRE(matrix.header.version == 2);
  REQUIRE(matrix.kmers.size() > 0);
  REQUIRE(matrix.patterns.size() > 0);
  REQUIRE(matrix.patterns.size() < matrix.kmers.size());

  // Spot-check multi-stripe bit packing: accession 65 uses stripe 1.
  bool found_stripe1 = false;
  for (const auto& words : matrix.patterns) {
    if (words.size() == 2 && words[1] != 0ULL) {
      found_stripe1 = true;
      break;
    }
  }
  REQUIRE(found_stripe1);

  kmat::ValidateOptions val_opts;
  val_opts.matrix_path = matrix_path;
  val_opts.accession_list_path = list_path;
  kmat::ValidateReport report;
  REQUIRE(kmat::validate_panel(val_opts, report).ok());
  REQUIRE(report.ok);

  const std::string pop_path = (out / "pop.tsv").string();
  kmat::PcaOptions pca_opts;
  pca_opts.matrix_path = matrix_path;
  pca_opts.accession_list_path = list_path;
  pca_opts.num_pcs = 0;  // all PCs
  kmat::PcaResult pca_result;
  REQUIRE(kmat::run_pca(pca_opts, pca_result).ok());
  REQUIRE(pca_result.accessions.size() == 72);
  REQUIRE(pca_result.pcs.front().size() == 72);
  REQUIRE(pca_result.accessions.front() == "acc_000");
  REQUIRE(kmat::write_pca_tsv(pop_path, pca_result).ok());

  kmat::GwasOptions gwas_opts;
  gwas_opts.matrix_path = matrix_path;
  gwas_opts.accession_list_path = list_path;
  gwas_opts.phenotype_path = (td / "phenotypes.tsv").string();
  gwas_opts.pop_path = pop_path;
  gwas_opts.kmer_size = 31;
  gwas_opts.num_pcs = 2;
  gwas_opts.print_all = true;
  kmat::GwasResult gwas_result;
  REQUIRE(kmat::run_gwas(gwas_opts, gwas_result).ok());
  REQUIRE(gwas_result.hits.size() > 0);

  kmat::GeneSearchOptions gene_opts;
  gene_opts.matrix_path = matrix_path;
  gene_opts.accession_list_path = list_path;
  gene_opts.gene_fasta_path = (td / "gene.fasta").string();
  gene_opts.kmer_size = 31;
  kmat::GeneSearchResult gene_result;
  REQUIRE(kmat::run_gene_search(gene_opts, gene_result).ok());
  REQUIRE(gene_result.hits.size() > 0);

  // Shared gene k-mers should be present in all 72 accessions.
  bool found_all_present = false;
  for (const auto& hit : gene_result.hits) {
    if (hit.pa_bits == std::string(72, '1')) {
      found_all_present = true;
      break;
    }
  }
  REQUIRE(found_all_present);
}
