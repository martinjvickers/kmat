#include "kmat/count.hpp"

#include "kmat/fastq.hpp"
#include "kmat/kmer.hpp"
#include "kmat/log.hpp"
#include "kmat/runtime.hpp"
#include "kmat/sequence.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

namespace kmat {

namespace {

Error validate_count_opts(const CountOptions& opts) {
  if (opts.input_path.empty() || opts.output_path.empty()) {
    return Error::invalid_argument("count requires input and output paths");
  }
  if (opts.kmer_size == 0 || opts.kmer_size > 32) {
    return Error::invalid_argument("k-mer size must be in 1..32");
  }
  if (opts.min_count == 0) {
    return Error::invalid_argument("min_count (--ci) must be >= 1");
  }
  return Error::success();
}

Error count_builtin(const CountOptions& opts) {
  std::unordered_map<std::uint64_t, std::uint32_t> counts;
  counts.reserve(1 << 20);

  if (auto err = for_each_sequence(opts.input_path, [&](const std::string& sequence) -> Error {
        return for_each_encoded_kmer(sequence, opts.kmer_size, [&](std::uint64_t code) {
          auto& c = counts[code];
          if (c < std::numeric_limits<std::uint32_t>::max()) {
            ++c;
          }
        });
      });
      !err.ok()) {
    return err;
  }

  PresenceSet set;
  set.header.kmer_size = static_cast<std::uint32_t>(opts.kmer_size);
  set.header.min_count = opts.min_count;
  set.kmers.reserve(counts.size());
  for (const auto& [code, c] : counts) {
    if (c >= opts.min_count) {
      set.kmers.push_back(code);
    }
  }
  std::sort(set.kmers.begin(), set.kmers.end());
  set.header.num_kmers = set.kmers.size();
  return write_presence_set(opts.output_path, set);
}

bool path_is_executable(const fs::path& p) {
  std::error_code ec;
  if (!fs::is_regular_file(p, ec)) {
    return false;
  }
  return ::access(p.c_str(), X_OK) == 0;
}

Error run_argv(const std::vector<std::string>& args, const std::string& stderr_path = {}) {
  if (args.empty()) {
    return Error::invalid_argument("empty command");
  }
  std::vector<char*> argv;
  argv.reserve(args.size() + 1);
  for (const std::string& a : args) {
    argv.push_back(const_cast<char*>(a.c_str()));
  }
  argv.push_back(nullptr);

  const pid_t pid = ::fork();
  if (pid < 0) {
    return Error::io_error("fork failed");
  }
  if (pid == 0) {
    if (!stderr_path.empty()) {
      const int fd = ::open(stderr_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
      if (fd >= 0) {
        ::dup2(fd, STDERR_FILENO);
        ::close(fd);
      }
    }
    ::execvp(argv[0], argv.data());
    ::_exit(127);
  }

  int status = 0;
  if (::waitpid(pid, &status, 0) < 0) {
    return Error::io_error("waitpid failed for " + args[0]);
  }
  if (WIFEXITED(status)) {
    const int code = WEXITSTATUS(status);
    if (code == 0) {
      return Error::success();
    }
    if (code == 127) {
      return Error::io_error("executable not found or not executable: " + args[0]);
    }
    std::ostringstream msg;
    msg << args[0] << " exited with status " << code;
    if (!stderr_path.empty()) {
      std::ifstream err_in(stderr_path);
      std::string err_line;
      std::string tail;
      while (std::getline(err_in, err_line)) {
        if (!err_line.empty()) {
          tail = err_line;
        }
      }
      if (!tail.empty()) {
        msg << " (" << tail << ")";
      }
    }
    return Error::io_error(msg.str());
  }
  return Error::io_error(args[0] + " terminated abnormally");
}

std::string scratch_root(const CountOptions& opts) {
  if (!opts.tmpdir.empty()) {
    return opts.tmpdir;
  }
  if (const char* t = std::getenv("TMPDIR")) {
    if (*t) {
      return t;
    }
  }
  return "/tmp";
}

const char* kmc_input_flag(const std::string& path) {
  if (path_looks_fasta(path)) {
    return "-fm";
  }
  // Default FASTQ (including .fq.gz — patched KMC in the image reads gzip via zlib gzFile).
  return "-fq";
}

Error dump_to_kset(const std::string& dump_path, const CountOptions& opts) {
  std::ifstream in(dump_path);
  if (!in) {
    return Error::io_error("failed to open KMC dump: " + dump_path);
  }

  PresenceSet set;
  set.header.kmer_size = static_cast<std::uint32_t>(opts.kmer_size);
  set.header.min_count = opts.min_count;
  set.kmers.reserve(1 << 20);

  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty()) {
      continue;
    }
    // Format: KMER<tab>COUNT  or  KMER COUNT
    std::string kmer;
    std::size_t i = 0;
    while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i]))) {
      kmer.push_back(line[i]);
      ++i;
    }
    if (kmer.size() != opts.kmer_size) {
      return Error::invalid_argument("KMC dump k-mer length mismatch: " + kmer);
    }
    std::uint64_t code = 0;
    if (auto err = encode_kmer(kmer, code); !err.ok()) {
      return err;
    }
    set.kmers.push_back(code);
  }

  std::sort(set.kmers.begin(), set.kmers.end());
  set.kmers.erase(std::unique(set.kmers.begin(), set.kmers.end()), set.kmers.end());
  set.header.num_kmers = set.kmers.size();
  return write_presence_set(opts.output_path, set);
}

Error count_kmc(const CountOptions& opts) {
  // Prefer patched KMC from the Singularity image (/usr/local/bin). Stock
  // release binaries use Cloudflare zlib + a manual inflate loop that fails on
  // many real Illumina/pigz .fq.gz files ("Some error while reading gzip file").
  std::string kmc_bin = opts.kmc_bin;
  std::string kmc_tools_bin = opts.kmc_tools_bin;
  if (kmc_bin.empty()) {
    if (path_is_executable("/usr/local/bin/kmc")) {
      kmc_bin = "/usr/local/bin/kmc";
    } else {
      kmc_bin = "kmc";
    }
  }
  if (kmc_tools_bin.empty()) {
    if (path_is_executable("/usr/local/bin/kmc_tools")) {
      kmc_tools_bin = "/usr/local/bin/kmc_tools";
    } else {
      kmc_tools_bin = "kmc_tools";
    }
  }

  std::string resolved_kmc;
  std::string resolved_tools;
  if (!find_executable(kmc_bin, resolved_kmc)) {
    return Error::io_error(
        "kmc not found on PATH (install patched KMC via singularity/kmat.def, or use "
        "--engine builtin). Looked for: " +
        kmc_bin);
  }
  if (!find_executable(kmc_tools_bin, resolved_tools)) {
    return Error::io_error(
        "kmc_tools not found on PATH (install patched KMC via singularity/kmat.def, or "
        "use --engine builtin). Looked for: " +
        kmc_tools_bin);
  }
  if (resolved_kmc == "/usr/bin/kmc") {
    log_info(
        "warning: using /usr/bin/kmc; prefer the patched KMC in /usr/local/bin from "
        "singularity/kmat.def (stock KMC often fails on .fq.gz)");
  }

  const std::size_t threads =
      opts.num_threads > 0 ? opts.num_threads : effective_threads(runtime_config());

  const fs::path root = fs::path(scratch_root(opts)) / ("kmat_kmc_" + std::to_string(::getpid()));
  std::error_code ec;
  fs::create_directories(root, ec);
  if (ec) {
    return Error::io_error("failed to create KMC work dir: " + root.string());
  }

  const fs::path db_prefix = root / "db";
  const fs::path work = root / "work";
  const fs::path dump = root / "dump.txt";
  const fs::path kmc_err = root / "kmc.err";
  const fs::path tools_err = root / "kmc_tools.err";
  fs::create_directories(work, ec);

  // RAII cleanup of scratch (KMC spill + dump only — never stages decompressed FASTQ).
  struct Cleaner {
    fs::path path;
    ~Cleaner() {
      std::error_code ignore;
      fs::remove_all(path, ignore);
    }
  } cleaner{root};

  const char* format_flag = kmc_input_flag(opts.input_path);

  std::ostringstream tflag;
  tflag << "-t" << threads;
  std::ostringstream kflag;
  kflag << "-k" << opts.kmer_size;
  std::ostringstream ciflag;
  ciflag << "-ci" << opts.min_count;
  std::ostringstream csflag;
  csflag << "-cs" << opts.min_count;

  log_info("KMC count: " + resolved_kmc + " " + kflag.str() + " " + ciflag.str() + " " +
           tflag.str() + " " + format_flag + " " + opts.input_path);

  {
    // Pass .fq.gz straight through. Image KMC is patched to gzread (no gunzip-to-disk).
    std::vector<std::string> args = {resolved_kmc,       "-hp",
                                     kflag.str(),        ciflag.str(),
                                     csflag.str(),       tflag.str(),
                                     format_flag,        opts.input_path,
                                     db_prefix.string(), work.string()};
    if (auto err = run_argv(args, kmc_err.string()); !err.ok()) {
      return err;
    }
  }

  {
    std::vector<std::string> args = {resolved_tools, "-hp", "transform", db_prefix.string(),
                                     "dump", "-s", dump.string()};
    if (auto err = run_argv(args, tools_err.string()); !err.ok()) {
      std::vector<std::string> alt = {resolved_tools, "-hp", "dump", "-s", db_prefix.string(),
                                      dump.string()};
      if (auto err2 = run_argv(alt, tools_err.string()); !err2.ok()) {
        return Error::io_error("kmc_tools dump failed (" + err.message + "; fallback: " +
                               err2.message + ")");
      }
    }
  }

  return dump_to_kset(dump.string(), opts);
}

}  // namespace

bool find_executable(const std::string& name, std::string& resolved_path) {
  resolved_path.clear();
  if (name.empty()) {
    return false;
  }
  const fs::path as_path(name);
  if (as_path.is_absolute() || name.find('/') != std::string::npos) {
    if (path_is_executable(as_path)) {
      resolved_path = as_path.string();
      return true;
    }
    return false;
  }

  const char* path_env = std::getenv("PATH");
  if (path_env == nullptr) {
    return false;
  }
  std::string path = path_env;
  std::size_t start = 0;
  while (start <= path.size()) {
    const std::size_t end = path.find(':', start);
    const std::string dir =
        path.substr(start, end == std::string::npos ? std::string::npos : end - start);
    if (!dir.empty()) {
      const fs::path candidate = fs::path(dir) / name;
      if (path_is_executable(candidate)) {
        resolved_path = candidate.string();
        return true;
      }
    }
    if (end == std::string::npos) {
      break;
    }
    start = end + 1;
  }
  return false;
}

Error count_kmers_to_presence_set(const CountOptions& opts) {
  if (auto err = validate_count_opts(opts); !err.ok()) {
    return err;
  }
  if (opts.engine == CountEngine::Builtin) {
    return count_builtin(opts);
  }
  return count_kmc(opts);
}

Error import_kmers_text_to_presence_set(const ImportKmersOptions& opts) {
  if (opts.input_path.empty() || opts.output_path.empty()) {
    return Error::invalid_argument("import-kmers requires input and output paths");
  }
  if (opts.kmer_size == 0 || opts.kmer_size > 32) {
    return Error::invalid_argument("k-mer size must be in 1..32");
  }

  std::ifstream in(opts.input_path);
  if (!in) {
    return Error::io_error("failed to open k-mer list: " + opts.input_path);
  }

  std::unordered_set<std::uint64_t> uniq;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty() || line.front() == '#') {
      continue;
    }
    while (!line.empty() && std::isspace(static_cast<unsigned char>(line.front()))) {
      line.erase(line.begin());
    }
    while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back()))) {
      line.pop_back();
    }
    if (line.empty()) {
      continue;
    }

    std::uint64_t code = 0;
    if (line.size() >= 2 && line[0] == '0' && (line[1] == 'x' || line[1] == 'X')) {
      try {
        code = std::stoull(line, nullptr, 16);
      } catch (...) {
        return Error::invalid_argument("invalid hex k-mer: " + line);
      }
    } else {
      if (line.size() != opts.kmer_size) {
        return Error::invalid_argument("k-mer length mismatch for: " + line);
      }
      if (auto err = encode_kmer(line, code); !err.ok()) {
        return err;
      }
    }
    uniq.insert(code);
  }

  PresenceSet set;
  set.header.kmer_size = static_cast<std::uint32_t>(opts.kmer_size);
  set.header.min_count = 1;
  set.kmers.assign(uniq.begin(), uniq.end());
  std::sort(set.kmers.begin(), set.kmers.end());
  set.header.num_kmers = set.kmers.size();
  return write_presence_set(opts.output_path, set);
}

}  // namespace kmat
