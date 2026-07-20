#include "kmat/count.hpp"
#include "kmat/gwas.hpp"
#include "kmat/io.hpp"
#include "kmat/matrix.hpp"
#include "kmat/pca.hpp"
#include "kmat/runtime.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Timing {
  std::string name;
  double ms{0};
};

double elapsed_ms(std::chrono::steady_clock::time_point t0,
                  std::chrono::steady_clock::time_point t1) {
  return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

void print_usage(const char* argv0) {
  std::cerr << "Usage: " << argv0
            << " --testdata DIR [--profile laptop|hpc] [--threads N] [--json OUT]\n"
            << "  Times build/pop/gwas on the medium panel (panel_k31_n72) or tiny panel.\n"
            << "  Use --panel tiny|medium (default medium).\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::string testdata;
  std::string profile_name = "laptop";
  std::size_t threads = 0;
  std::string json_out;
  std::string panel = "medium";

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto need = [&](const char* flag) -> std::string {
      if (i + 1 >= argc) {
        std::cerr << "missing value for " << flag << '\n';
        std::exit(2);
      }
      return argv[++i];
    };
    if (arg == "--testdata") {
      testdata = need("--testdata");
    } else if (arg == "--profile") {
      profile_name = need("--profile");
    } else if (arg == "--threads") {
      threads = static_cast<std::size_t>(std::stoul(need("--threads")));
    } else if (arg == "--json") {
      json_out = need("--json");
    } else if (arg == "--panel") {
      panel = need("--panel");
    } else if (arg == "-h" || arg == "--help") {
      print_usage(argv[0]);
      return 0;
    } else {
      std::cerr << "unknown arg: " << arg << '\n';
      print_usage(argv[0]);
      return 2;
    }
  }

  if (testdata.empty()) {
    print_usage(argv[0]);
    return 2;
  }

  bool ok = false;
  const kmat::RuntimeProfile profile = kmat::parse_runtime_profile(profile_name, ok);
  if (!ok) {
    std::cerr << "unknown profile: " << profile_name << '\n';
    return 2;
  }
  const kmat::RuntimeConfig cfg = kmat::resolve_runtime(profile, threads);
  kmat::set_runtime_config(cfg);

  fs::path td = fs::path(testdata);
  std::size_t k = 31;
  fs::path accession_list;
  fs::path phenotypes;
  if (panel == "medium") {
    if (fs::is_directory(td / "panel_k31_n72")) {
      td = td / "panel_k31_n72";
    }
    accession_list = td / "accession_list.txt";
    phenotypes = td / "phenotypes.tsv";
    k = 31;
  } else if (panel == "tiny") {
    accession_list = td / "accession_list.txt";
    phenotypes = td / "phenotypes.tsv";
    k = 3;
  } else {
    std::cerr << "unknown --panel (use tiny or medium)\n";
    return 2;
  }

  if (!fs::exists(accession_list) || !fs::exists(phenotypes)) {
    std::cerr << "testdata paths missing under " << td << '\n';
    return 2;
  }

  const fs::path work = fs::temp_directory_path() / "kmat_bench";
  fs::create_directories(work);
  const fs::path matrix_path = work / (panel + "_" + profile_name + ".kmat");
  const fs::path pop_path = work / (panel + "_" + profile_name + "_pop.tsv");

  std::vector<std::string> paths;
  if (auto err = kmat::read_list_file(accession_list.string(), paths); !err.ok()) {
    std::cerr << err.message << '\n';
    return 1;
  }
  if (auto err = kmat::resolve_list_paths(accession_list.string(), paths); !err.ok()) {
    std::cerr << err.message << '\n';
    return 1;
  }

  std::vector<Timing> timings;
  const auto threads_eff = kmat::effective_threads(cfg);

  {
    kmat::BuildOptions opts;
    opts.kmer_size = k;
    opts.accession_paths = paths;
    opts.output_path = matrix_path.string();
    opts.num_threads = threads_eff;
    const auto t0 = std::chrono::steady_clock::now();
    if (auto err = kmat::build_matrix_from_accessions(opts); !err.ok()) {
      std::cerr << "build failed: " << err.message << '\n';
      return 1;
    }
    timings.push_back({"build", elapsed_ms(t0, std::chrono::steady_clock::now())});
  }

  {
    kmat::PcaOptions opts;
    opts.matrix_path = matrix_path.string();
    opts.accession_list_path = accession_list.string();
    opts.num_pcs = 2;
    const auto t0 = std::chrono::steady_clock::now();
    kmat::PcaResult result;
    if (auto err = kmat::run_pca(opts, result); !err.ok()) {
      std::cerr << "pop failed: " << err.message << '\n';
      return 1;
    }
    if (auto err = kmat::write_pca_tsv(pop_path.string(), result); !err.ok()) {
      std::cerr << "pop write failed: " << err.message << '\n';
      return 1;
    }
    timings.push_back({"pop", elapsed_ms(t0, std::chrono::steady_clock::now())});
  }

  {
    kmat::GwasOptions opts;
    opts.matrix_path = matrix_path.string();
    opts.accession_list_path = accession_list.string();
    opts.phenotype_path = phenotypes.string();
    opts.pop_path = pop_path.string();
    opts.kmer_size = k;
    opts.print_all = true;
    opts.num_threads = threads_eff;
    const auto t0 = std::chrono::steady_clock::now();
    kmat::GwasResult result;
    if (auto err = kmat::run_gwas(opts, result); !err.ok()) {
      std::cerr << "gwas failed: " << err.message << '\n';
      return 1;
    }
    timings.push_back({"gwas", elapsed_ms(t0, std::chrono::steady_clock::now())});
  }

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "kmat_bench profile=" << kmat::runtime_profile_name(profile)
            << " threads=" << threads_eff << " panel=" << panel << '\n';
  for (const Timing& t : timings) {
    std::cout << "  " << t.name << "_ms=" << t.ms << '\n';
  }

  if (!json_out.empty()) {
    std::ofstream out(json_out);
    if (!out) {
      std::cerr << "failed to write " << json_out << '\n';
      return 1;
    }
    out << "{\n";
    out << "  \"profile\": \"" << kmat::runtime_profile_name(profile) << "\",\n";
    out << "  \"threads\": " << threads_eff << ",\n";
    out << "  \"panel\": \"" << panel << "\",\n";
    out << "  \"io_buffer_bytes\": " << cfg.io_buffer_bytes << ",\n";
    out << "  \"timings_ms\": {\n";
    for (std::size_t i = 0; i < timings.size(); ++i) {
      out << "    \"" << timings[i].name << "\": " << timings[i].ms;
      if (i + 1 < timings.size()) {
        out << ',';
      }
      out << '\n';
    }
    out << "  }\n";
    out << "}\n";
  }

  return 0;
}
