#include "kmat/count.hpp"
#include "kmat/gene.hpp"
#include "kmat/gwas.hpp"
#include "kmat/io.hpp"
#include "kmat/log.hpp"
#include "kmat/matrix.hpp"
#include "kmat/pca.hpp"
#include "kmat/runtime.hpp"
#include "kmat/stripe_build.hpp"
#include "kmat/universe.hpp"
#include "kmat/validate.hpp"

#include <CLI/CLI.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
namespace {

int fail(const kmat::Error& err) {
  kmat::log_error(err.message);
  std::cerr << "error: " << err.message << '\n';
  return 1;
}

int run_count(const std::string& input, std::size_t kmer_size, std::uint32_t min_count,
              const std::string& output, const std::string& engine_name, const std::string& tmpdir) {
  kmat::CountOptions opts;
  opts.input_path = input;
  opts.output_path = output;
  opts.kmer_size = kmer_size;
  opts.min_count = min_count;
  opts.tmpdir = tmpdir;
  opts.num_threads = 0;  // use runtime profile
  if (engine_name == "builtin") {
    opts.engine = kmat::CountEngine::Builtin;
  } else if (engine_name == "kmc") {
    opts.engine = kmat::CountEngine::Kmc;
  } else {
    std::cerr << "error: --engine must be kmc or builtin\n";
    return 1;
  }
  if (auto err = kmat::count_kmers_to_presence_set(opts); !err.ok()) {
    return fail(err);
  }
  kmat::log_info("wrote presence set: " + output);
  return 0;
}

int run_import_kmers(const std::string& input, std::size_t kmer_size, const std::string& output) {
  kmat::ImportKmersOptions opts;
  opts.input_path = input;
  opts.output_path = output;
  opts.kmer_size = kmer_size;
  if (auto err = kmat::import_kmers_text_to_presence_set(opts); !err.ok()) {
    return fail(err);
  }
  kmat::log_info("wrote presence set: " + output);
  return 0;
}

int run_build(const std::string& accession_list, std::size_t kmer_size, const std::string& output,
              double memory_gb, std::size_t batch_rows, const std::string& tmpdir) {
  std::vector<std::string> paths;
  if (auto err = kmat::read_list_file(accession_list, paths); !err.ok()) {
    return fail(err);
  }
  if (auto err = kmat::resolve_list_paths(accession_list, paths); !err.ok()) {
    return fail(err);
  }

  kmat::BuildOptions opts;
  opts.kmer_size = kmer_size;
  opts.accession_paths = std::move(paths);
  opts.output_path = output;
  opts.batch_rows = batch_rows;
  opts.tmpdir = tmpdir;
  if (memory_gb > 0.0) {
    opts.memory_bytes = static_cast<std::size_t>(memory_gb * (1024.0 * 1024.0 * 1024.0));
  }

  if (auto err = kmat::build_matrix_from_accessions(opts); !err.ok()) {
    return fail(err);
  }

  kmat::log_info("wrote PA matrix: " + output);
  return 0;
}

int run_build_master(const std::string& accession_list, std::size_t kmer_size,
                     const std::string& output, const std::string& tmpdir, std::size_t group_size) {
  std::vector<std::string> paths;
  if (auto err = kmat::read_list_file(accession_list, paths); !err.ok()) {
    return fail(err);
  }
  if (auto err = kmat::resolve_list_paths(accession_list, paths); !err.ok()) {
    return fail(err);
  }
  if (auto err = kmat::build_universe_from_presence_sets(paths, kmer_size, output, tmpdir,
                                                         group_size > 0 ? group_size : 32);
      !err.ok()) {
    return fail(err);
  }
  kmat::log_info("wrote universe: " + output);
  return 0;
}

int run_create_stripe(const std::string& universe, const std::string& output, std::size_t kmer_size,
                      std::size_t num_accessions, std::size_t stripe_index, std::size_t batch_rows) {
  kmat::CreateStripeOptions opts;
  opts.universe_path = universe;
  opts.output_path = output;
  opts.kmer_size = kmer_size;
  opts.num_accessions = num_accessions;
  opts.stripe_index = stripe_index;
  opts.batch_rows = batch_rows;
  if (auto err = kmat::create_blank_stripe(opts); !err.ok()) {
    return fail(err);
  }
  kmat::log_info("wrote blank stripe: " + output);
  return 0;
}

int run_fill(const std::string& stripe, const std::string& kset, std::size_t local_index,
             std::size_t batch_rows) {
  kmat::FillStripeOptions opts;
  opts.stripe_path = stripe;
  opts.kset_path = kset;
  opts.local_index = local_index;
  opts.batch_rows = batch_rows;
  if (auto err = kmat::fill_stripe_column(opts); !err.ok()) {
    return fail(err);
  }
  kmat::log_info("filled stripe column: " + stripe + " n=" + std::to_string(local_index));
  return 0;
}

int run_compress(const std::string& matrix_list, const std::string& accession_list,
                 std::size_t kmer_size, const std::string& output, double memory_gb,
                 std::size_t batch_rows, const std::string& tmpdir) {
  std::vector<std::string> stripes;
  if (auto err = kmat::read_matrix_list(matrix_list, stripes); !err.ok()) {
    return fail(err);
  }
  std::vector<std::string> accessions;
  if (auto err = kmat::read_list_file(accession_list, accessions); !err.ok()) {
    return fail(err);
  }
  kmat::CompressOptions opts;
  opts.stripe_paths = std::move(stripes);
  opts.num_accessions = accessions.size();
  opts.kmer_size = kmer_size;
  opts.output_path = output;
  opts.batch_rows = batch_rows;
  opts.tmpdir = tmpdir;
  if (memory_gb > 0.0) {
    opts.memory_bytes = static_cast<std::size_t>(memory_gb * (1024.0 * 1024.0 * 1024.0));
  }
  if (auto err = kmat::compress_stripes_to_v2(opts); !err.ok()) {
    return fail(err);
  }
  kmat::log_info("wrote v2 matrix: " + output);
  return 0;
}

int run_pop(const std::string& matrix_path, const std::string& matrix_list,
            const std::string& accession_list, const std::string& output, std::size_t num_pcs,
            std::size_t max_samples, unsigned seed) {
  kmat::PcaOptions opts;
  opts.matrix_path = matrix_path;
  opts.matrix_list_path = matrix_list;
  opts.accession_list_path = accession_list;
  opts.num_pcs = num_pcs;
  opts.max_samples = max_samples;
  opts.seed = seed;

  kmat::PcaResult result;
  if (auto err = kmat::run_pca(opts, result); !err.ok()) {
    return fail(err);
  }
  if (auto err = kmat::write_pca_tsv(output, result); !err.ok()) {
    return fail(err);
  }

  kmat::log_info("wrote population structure: " + output);
  return 0;
}

int run_gwas(const std::string& matrix_path, const std::string& matrix_list,
             const std::string& accession_list, const std::string& phenotype_path,
             const std::string& pop_path, std::size_t kmer_size, std::size_t num_pcs,
             std::size_t top_n, bool print_all, bool display_pa) {
  kmat::GwasOptions opts;
  opts.matrix_path = matrix_path;
  opts.matrix_list_path = matrix_list;
  opts.accession_list_path = accession_list;
  opts.phenotype_path = phenotype_path;
  opts.pop_path = pop_path;
  opts.kmer_size = kmer_size;
  opts.num_pcs = num_pcs;
  opts.top_n = top_n;
  opts.print_all = print_all;
  opts.include_pa_bits = display_pa;

  kmat::GwasResult result;
  if (auto err = kmat::run_gwas(opts, result); !err.ok()) {
    return fail(err);
  }
  if (auto err = kmat::write_gwas_tsv(std::cout, kmer_size, result, display_pa); !err.ok()) {
    return fail(err);
  }
  return 0;
}

int run_gene(const std::string& matrix_path, const std::string& matrix_list,
             const std::string& accession_list, const std::string& gene_fasta, std::size_t kmer_size) {
  kmat::GeneSearchOptions opts;
  opts.matrix_path = matrix_path;
  opts.matrix_list_path = matrix_list;
  opts.accession_list_path = accession_list;
  opts.gene_fasta_path = gene_fasta;
  opts.kmer_size = kmer_size;

  kmat::GeneSearchResult result;
  if (auto err = kmat::run_gene_search(opts, result); !err.ok()) {
    return fail(err);
  }
  if (auto err = kmat::write_gene_hits(std::cout, result); !err.ok()) {
    return fail(err);
  }
  return 0;
}

int run_validate(const std::string& matrix_path, const std::string& matrix_list,
                 const std::string& accession_list) {
  kmat::ValidateOptions opts;
  opts.matrix_path = matrix_path;
  opts.matrix_list_path = matrix_list;
  opts.accession_list_path = accession_list;

  kmat::ValidateReport report;
  if (auto err = kmat::validate_panel(opts, report); !err.ok()) {
    return fail(err);
  }

  if (report.ok) {
    std::cout << report.summary << '\n';
    return 0;
  }

  std::cerr << report.summary << '\n';
  for (const std::string& issue : report.issues) {
    std::cerr << "  - " << issue << '\n';
  }
  return 1;
}

}  // namespace

int main(int argc, char** argv) {
  CLI::App app{"kmat — k-mer presence/absence GWAS toolkit"};
  app.require_subcommand(0, 1);
  app.set_version_flag("-V,--version", "0.4.0");
  bool verbose = false;
  app.add_flag("-v,--verbose", verbose, "Enable verbose logging");

  std::string profile_name = "laptop";
  std::size_t threads_override = 0;
  app.add_option("--profile", profile_name, "Runtime profile: laptop|hpc")->default_val("laptop");
  app.add_option("--threads", threads_override, "Override worker thread count (0=profile default)")
      ->default_val(0);

  std::string count_input;
  std::size_t count_k = 31;
  std::uint32_t count_ci = 1;
  std::string count_output;
  std::string count_engine = "kmc";
  std::string count_tmpdir;
  CLI::App* count = app.add_subcommand("count", "Count/filter k-mers from FASTQ/FASTA into a .kset");
  count->add_option("-i,--input", count_input, "Input FASTQ/FASTA (optional .gz)")->required();
  count->add_option("-s,--kmer-size", count_k, "K-mer size")->required();
  count->add_option("--ci", count_ci, "Minimum k-mer count (KMC-style -ci)")->default_val(1);
  count->add_option("-o,--output", count_output, "Output .kset presence set")->required();
  count->add_option("--engine", count_engine, "Count engine: kmc (default) or builtin")
      ->default_val("kmc");
  count->add_option("--tmpdir", count_tmpdir, "Scratch directory for KMC (default: $TMPDIR)");

  std::string import_input;
  std::size_t import_k = 31;
  std::string import_output;
  CLI::App* import_kmers =
      app.add_subcommand("import-kmers", "Import a text k-mer list into a .kset (KMC migration)");
  import_kmers->add_option("-i,--input", import_input, "Text file: one DNA k-mer per line")->required();
  import_kmers->add_option("-s,--kmer-size", import_k, "K-mer size")->required();
  import_kmers->add_option("-o,--output", import_output, "Output .kset")->required();

  std::string build_accession_list;
  std::size_t build_k = 31;
  std::string build_output;
  double build_memory_gb = 0.0;
  std::size_t build_batch_rows = 0;
  std::string build_tmpdir;
  CLI::App* build =
      app.add_subcommand("build", "Build a PA matrix from sequences or .kset presence sets");
  build->add_option("-k,--accession-list", build_accession_list,
                    "Accession list (FASTA/FASTQ/.gz or .kset paths)")
      ->required();
  build->add_option("-s,--kmer-size", build_k, "K-mer size")->required();
  build->add_option("-o,--output", build_output, "Output matrix path")->required();
  build->add_option("--memory-gb", build_memory_gb,
                    "Working-set budget in GiB for streaming .kset build (0=profile default)");
  build->add_option("--batch-rows", build_batch_rows,
                    "Rows per create/fill I/O batch (0=100000)");
  build->add_option("--tmpdir", build_tmpdir, "Scratch directory for staged build");

  std::string master_list;
  std::size_t master_k = 31;
  std::string master_output;
  std::string master_tmpdir;
  std::size_t master_group = 32;
  CLI::App* build_master =
      app.add_subcommand("build-master", "Tree-merge .kset files into a .kuniv master universe");
  build_master->add_option("-k,--accession-list", master_list, "List of .kset paths")->required();
  build_master->add_option("-s,--kmer-size", master_k, "K-mer size")->required();
  build_master->add_option("-o,--output", master_output, "Output .kuniv path")->required();
  build_master->add_option("--tmpdir", master_tmpdir, "Scratch for merge temps");
  build_master->add_option("--group-size", master_group, "Tree-merge fan-in")->default_val(32);

  std::string create_universe;
  std::string create_output;
  std::size_t create_k = 31;
  std::size_t create_nacc = 0;
  std::size_t create_stripe = 0;
  std::size_t create_batch = 0;
  CLI::App* create_stripe_cmd =
      app.add_subcommand("create-stripe", "Create a blank v1 dense stripe from .kuniv");
  create_stripe_cmd->add_option("-m,--universe", create_universe, "Master .kuniv")->required();
  create_stripe_cmd->add_option("-o,--output", create_output, "Output stripe .bin")->required();
  create_stripe_cmd->add_option("-s,--kmer-size", create_k, "K-mer size")->required();
  create_stripe_cmd->add_option("--num-accessions", create_nacc, "Global accession count")
      ->required();
  create_stripe_cmd->add_option("--stripe-index", create_stripe, "Stripe index (0-based)")
      ->required();
  create_stripe_cmd->add_option("--batch-rows", create_batch, "I/O batch (0=100000)");

  std::string fill_stripe;
  std::string fill_kset;
  std::string fill_list;
  std::size_t fill_n = 0;
  std::size_t fill_batch = 0;
  CLI::App* fill = app.add_subcommand("fill", "Fill one local column in a dense stripe from a .kset");
  fill->add_option("-o,--output", fill_stripe, "Stripe .bin to update")->required();
  fill->add_option("--kset", fill_kset, "Accession .kset path");
  fill->add_option("-k,--accession-list", fill_list, "Stripe-local list (use with -n)");
  fill->add_option("-n,--local-index", fill_n, "Local bit index within stripe (0..63)")->required();
  fill->add_option("--batch-rows", fill_batch, "I/O batch (0=100000)");

  std::string compress_matrix_list;
  std::string compress_accession_list;
  std::size_t compress_k = 31;
  std::string compress_output;
  double compress_memory_gb = 0.0;
  std::size_t compress_batch = 0;
  std::string compress_tmpdir;
  CLI::App* compress =
      app.add_subcommand("compress", "Compress filled v1 stripes into a v2 .kmat");
  compress->add_option("-m,--matrix-list", compress_matrix_list, "List of filled stripe .bin paths")
      ->required();
  compress->add_option("-k,--accession-list", compress_accession_list, "Global accession list")
      ->required();
  compress->add_option("-s,--kmer-size", compress_k, "K-mer size")->required();
  compress->add_option("-o,--output", compress_output, "Output v2 .kmat")->required();
  compress->add_option("--memory-gb", compress_memory_gb, "Pattern-dict budget (0=profile default)");
  compress->add_option("--batch-rows", compress_batch, "I/O batch (0=100000)");
  compress->add_option("--tmpdir", compress_tmpdir, "Scratch directory");

  std::string pop_matrix;
  std::string pop_matrix_list;
  std::string pop_accession_list;
  std::string pop_output;
  std::size_t pop_npc = 0;
  std::size_t pop_max_samples = 0;
  unsigned pop_seed = 42;
  CLI::App* pop = app.add_subcommand("pop", "Compute population-structure covariates (PCA)");
  pop->add_option("-i,--matrix", pop_matrix, "Single PA matrix file");
  pop->add_option("-m,--matrix-list", pop_matrix_list, "PA stripe list (matrix_list.txt)");
  pop->add_option("-k,--accession-list", pop_accession_list, "Accession list")->required();
  pop->add_option("-o,--output", pop_output, "Output pop TSV")->required();
  pop->add_option("--npc", pop_npc, "Number of PCs to emit (0=all)")->default_val(0);
  pop->add_option("--max-samples", pop_max_samples,
                  "Max patterns to sample for PCA (0=auto: all if ≤500k else 100k)");
  pop->add_option("--seed", pop_seed, "RNG seed for sampling");

  std::string gwas_matrix;
  std::string gwas_matrix_list;
  std::string gwas_accession_list;
  std::string gwas_phenotype;
  std::string gwas_pop;
  std::size_t gwas_k = 31;
  std::size_t gwas_npc = 2;
  std::size_t gwas_top_n = 1000;
  bool gwas_print_all = false;
  bool gwas_display = false;
  CLI::App* gwas = app.add_subcommand("gwas", "Run genome-wide k-mer association");
  gwas->add_option("-i,--matrix", gwas_matrix, "Single PA matrix file");
  gwas->add_option("-m,--matrix-list", gwas_matrix_list, "PA stripe list");
  gwas->add_option("-k,--accession-list", gwas_accession_list, "Accession list")->required();
  gwas->add_option("-p,--phenotype-file", gwas_phenotype, "Phenotype TSV")->required();
  gwas->add_option("--pop,--population-structure-file", gwas_pop, "Pop structure TSV")->required();
  gwas->add_option("-s,--kmer-size", gwas_k, "K-mer size")->required();
  gwas->add_option("--npc", gwas_npc, "Number of leading PCs from pop TSV to use as covariates")
      ->default_val(2);
  gwas->add_option("-N,--top-n", gwas_top_n, "Keep top N hits by min p-value")->default_val(1000);
  gwas->add_flag("-a,--print-all", gwas_print_all, "Score and print all k-mers");
  gwas->add_flag("-d,--display", gwas_display, "Include PA bitstring column");

  std::string gene_matrix;
  std::string gene_matrix_list;
  std::string gene_accession_list;
  std::string gene_fasta;
  std::size_t gene_k = 31;
  CLI::App* gene = app.add_subcommand("gene", "Look up PA patterns for gene k-mers");
  gene->add_option("-i,--matrix", gene_matrix, "Single PA matrix file");
  gene->add_option("-m,--matrix-list", gene_matrix_list, "PA stripe list");
  gene->add_option("-k,--accession-list", gene_accession_list, "Accession list")->required();
  gene->add_option("-g,--gene", gene_fasta, "Gene FASTA")->required();
  gene->add_option("-s,--kmer-size", gene_k, "K-mer size")->required();

  std::string validate_matrix;
  std::string validate_matrix_list;
  std::string validate_accession_list;
  CLI::App* validate = app.add_subcommand("validate", "Validate matrix and list-file consistency");
  validate->add_option("-i,--matrix", validate_matrix, "Single PA matrix file");
  validate->add_option("-m,--matrix-list", validate_matrix_list, "PA stripe list");
  validate->add_option("-k,--accession-list", validate_accession_list, "Accession list")->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError& e) {
    return app.exit(e);
  }

  if (verbose) {
    kmat::set_log_level(kmat::LogLevel::Debug);
  }

  if (app.get_subcommands().empty()) {
    std::cout << app.help() << '\n';
    return 0;
  }

  bool profile_ok = false;
  const kmat::RuntimeProfile profile = kmat::parse_runtime_profile(profile_name, profile_ok);
  if (!profile_ok) {
    std::cerr << "error: unknown --profile '" << profile_name << "' (use laptop or hpc)\n";
    return 1;
  }
  kmat::set_runtime_config(kmat::resolve_runtime(profile, threads_override));
  kmat::log_info(std::string("runtime profile=") + kmat::runtime_profile_name(profile) +
                 " threads=" + std::to_string(kmat::effective_threads(kmat::runtime_config())));

  CLI::App* sub = app.get_subcommands().front();
  const std::string name = sub->get_name();

  if (name == "count") {
    return run_count(count_input, count_k, count_ci, count_output, count_engine, count_tmpdir);
  }
  if (name == "import-kmers") {
    return run_import_kmers(import_input, import_k, import_output);
  }
  if (name == "build") {
    return run_build(build_accession_list, build_k, build_output, build_memory_gb, build_batch_rows,
                     build_tmpdir);
  }
  if (name == "build-master") {
    return run_build_master(master_list, master_k, master_output, master_tmpdir, master_group);
  }
  if (name == "create-stripe") {
    return run_create_stripe(create_universe, create_output, create_k, create_nacc, create_stripe,
                             create_batch);
  }
  if (name == "fill") {
    std::string kset = fill_kset;
    if (kset.empty()) {
      if (fill_list.empty()) {
        std::cerr << "error: provide --kset or -k accession-list with -n\n";
        return 1;
      }
      std::vector<std::string> paths;
      if (auto err = kmat::read_list_file(fill_list, paths); !err.ok()) {
        return fail(err);
      }
      if (auto err = kmat::resolve_list_paths(fill_list, paths); !err.ok()) {
        return fail(err);
      }
      if (fill_n >= paths.size()) {
        std::cerr << "error: -n out of range for list\n";
        return 1;
      }
      kset = paths[fill_n];
    }
    return run_fill(fill_stripe, kset, fill_n, fill_batch);
  }
  if (name == "compress") {
    return run_compress(compress_matrix_list, compress_accession_list, compress_k, compress_output,
                        compress_memory_gb, compress_batch, compress_tmpdir);
  }
  if (name == "pop") {
    if (pop_matrix.empty() && pop_matrix_list.empty()) {
      std::cerr << "error: provide --matrix or --matrix-list\n";
      return 1;
    }
    return run_pop(pop_matrix, pop_matrix_list, pop_accession_list, pop_output, pop_npc,
                   pop_max_samples, pop_seed);
  }
  if (name == "gwas") {
    if (gwas_matrix.empty() && gwas_matrix_list.empty()) {
      std::cerr << "error: provide --matrix or --matrix-list\n";
      return 1;
    }
    return run_gwas(gwas_matrix, gwas_matrix_list, gwas_accession_list, gwas_phenotype, gwas_pop,
                    gwas_k, gwas_npc, gwas_top_n, gwas_print_all, gwas_display);
  }
  if (name == "gene") {
    if (gene_matrix.empty() && gene_matrix_list.empty()) {
      std::cerr << "error: provide --matrix or --matrix-list\n";
      return 1;
    }
    return run_gene(gene_matrix, gene_matrix_list, gene_accession_list, gene_fasta, gene_k);
  }
  if (name == "validate") {
    if (validate_matrix.empty() && validate_matrix_list.empty()) {
      std::cerr << "error: provide --matrix or --matrix-list\n";
      return 1;
    }
    return run_validate(validate_matrix, validate_matrix_list, validate_accession_list);
  }

  std::cerr << "kmat: unknown subcommand\n";
  return 1;
}
