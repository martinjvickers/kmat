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

fs::path temp_dir() {
  const fs::path base = fs::temp_directory_path() / "kmat_integration";
  std::error_code ec;
  fs::create_directories(base, ec);
  return base;
}

}  // namespace

TEST_CASE("end-to-end build pop gwas gene on testdata", "[integration]") {
  const fs::path td = fs::path(KMAT_TESTDATA_DIR);
  const fs::path out = temp_dir();

  std::vector<std::string> paths;
  const std::string list_path = (td / "accession_list.txt").string();
  REQUIRE(kmat::read_list_file(list_path, paths).ok());
  REQUIRE(kmat::resolve_list_paths(list_path, paths).ok());

  const std::string matrix_path = (out / "panel.kmat").string();
  kmat::BuildOptions build_opts;
  build_opts.kmer_size = 3;
  build_opts.accession_paths = paths;
  build_opts.output_path = matrix_path;
  REQUIRE(kmat::build_matrix_from_sequences(build_opts).ok());

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
  pca_opts.num_pcs = 2;
  kmat::PcaResult pca_result;
  REQUIRE(kmat::run_pca(pca_opts, pca_result).ok());
  REQUIRE(pca_result.accessions.size() == 6);
  REQUIRE(pca_result.pcs.size() == 6);
  REQUIRE(pca_result.pcs.front().size() == 2);
  REQUIRE(kmat::write_pca_tsv(pop_path, pca_result).ok());

  kmat::GwasOptions gwas_opts;
  gwas_opts.matrix_path = matrix_path;
  gwas_opts.accession_list_path = list_path;
  gwas_opts.phenotype_path = (td / "phenotypes.tsv").string();
  gwas_opts.pop_path = pop_path;
  gwas_opts.kmer_size = 3;
  gwas_opts.print_all = true;
  kmat::GwasResult gwas_result;
  REQUIRE(kmat::run_gwas(gwas_opts, gwas_result).ok());
  REQUIRE(gwas_result.hits.size() > 0);

  std::ostringstream gwas_out;
  REQUIRE(kmat::write_gwas_tsv(gwas_out, 3, gwas_result, false).ok());
  REQUIRE(gwas_out.str().find("kmer\t") != std::string::npos);

  kmat::GeneSearchOptions gene_opts;
  gene_opts.matrix_path = matrix_path;
  gene_opts.accession_list_path = list_path;
  gene_opts.gene_fasta_path = (td / "gene.fasta").string();
  gene_opts.kmer_size = 3;
  kmat::GeneSearchResult gene_result;
  REQUIRE(kmat::run_gene_search(gene_opts, gene_result).ok());
  REQUIRE(gene_result.hits.size() > 0);

  std::ostringstream gene_out;
  REQUIRE(kmat::write_gene_hits(gene_out, gene_result).ok());
  REQUIRE(gene_out.str().find("ACG") != std::string::npos);
}
