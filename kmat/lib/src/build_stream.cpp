#include "kmat/build_stream.hpp"

#include "kmat/log.hpp"
#include "kmat/matrix_layout.hpp"
#include "kmat/presence.hpp"
#include "kmat/runtime.hpp"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <numeric>
#include <optional>
#include <queue>
#include <sstream>
#include <string>
#include <sys/resource.h>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace kmat {

namespace {

struct PatternKey {
  std::vector<std::uint64_t> words;
  bool operator==(const PatternKey& o) const { return words == o.words; }
};

struct PatternKeyHash {
  std::size_t operator()(const PatternKey& k) const noexcept {
    std::size_t h = k.words.size();
    for (std::uint64_t w : k.words) {
      h ^= std::hash<std::uint64_t>{}(w) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    }
    return h;
  }
};

std::size_t code_partition(std::uint64_t code, std::size_t num_partitions) {
  std::uint64_t x = code + 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  x = x ^ (x >> 31);
  return static_cast<std::size_t>(x % num_partitions);
}

std::string scratch_root(const BuildOptions& opts) {
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

std::size_t default_memory_bytes() {
  const auto cfg = runtime_config();
  if (cfg.profile == RuntimeProfile::Hpc) {
    return 64ull << 30;
  }
  return 8ull << 30;
}

void raise_nofile_soft_to_hard() {
  struct rlimit rl {};
  if (::getrlimit(RLIMIT_NOFILE, &rl) != 0) {
    return;
  }
  if (rl.rlim_cur < rl.rlim_max) {
    rl.rlim_cur = rl.rlim_max;
    (void)::setrlimit(RLIMIT_NOFILE, &rl);
  }
}

std::size_t merge_fd_budget() {
  if (const char* e = std::getenv("KMAT_BUILD_MAX_OPEN")) {
    char* end = nullptr;
    const unsigned long v = std::strtoul(e, &end, 10);
    if (end != e && v >= 2) {
      return static_cast<std::size_t>(v);
    }
  }
  raise_nofile_soft_to_hard();
  std::size_t nofile = 1024;
  struct rlimit rl {};
  if (::getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur > 0 &&
      rl.rlim_cur != RLIM_INFINITY) {
    nofile = static_cast<std::size_t>(rl.rlim_cur);
  }
  constexpr std::size_t kReserve = 64;
  if (nofile <= kReserve + 2) {
    return 2;
  }
  return nofile - kReserve;
}

std::size_t choose_partitions(std::size_t n_stripes, std::size_t memory_bytes,
                              std::size_t num_threads, std::size_t batch_rows,
                              std::uint64_t sum_kmers) {
  const std::size_t threads = std::max<std::size_t>(1, num_threads);
  const std::size_t record_bytes = sizeof(std::uint64_t) * (1 + n_stripes);
  const std::size_t batch = std::max<std::size_t>(1, batch_rows);

  const std::size_t budget_per_thread =
      std::max<std::size_t>(1ull << 20, memory_bytes / threads);
  const std::uint64_t u_upper = std::max<std::uint64_t>(1, sum_kmers);
  const std::uint64_t bytes_per_row_est = record_bytes * 2;
  const std::uint64_t rows_per_shard =
      std::max<std::uint64_t>(1, budget_per_thread / bytes_per_row_est);
  std::uint64_t p64 = (u_upper + rows_per_shard - 1) / rows_per_shard;
  p64 = std::max<std::uint64_t>(threads, p64);

  const std::size_t buf_budget = std::max<std::size_t>(record_bytes * batch, memory_bytes / 4);
  const std::uint64_t max_p_by_buf =
      std::max<std::uint64_t>(1, buf_budget / (record_bytes * batch));
  p64 = std::min(p64, max_p_by_buf);

  constexpr std::uint64_t kAbsMaxPartitions = 1ull << 16;
  p64 = std::min(p64, kAbsMaxPartitions);
  p64 = std::max<std::uint64_t>(1, p64);
  return static_cast<std::size_t>(p64);
}

std::size_t scatter_batch_rows(std::size_t batch_rows, std::size_t num_partitions,
                               std::size_t num_threads, std::size_t memory_bytes) {
  const std::size_t batch = std::max<std::size_t>(1, batch_rows);
  const std::size_t p = std::max<std::size_t>(1, num_partitions);
  const std::size_t t = std::max<std::size_t>(1, num_threads);
  // One worker may buffer T code buckets; keep aggregate under ~1/8 memory.
  const std::size_t budget = std::max<std::size_t>(8, memory_bytes / 8);
  const std::size_t max_batch = std::max<std::size_t>(1, budget / (t * p * sizeof(std::uint64_t)));
  return std::min(batch, max_batch);
}

int progress_log_interval_sec() {
  if (const char* e = std::getenv("KMAT_BUILD_LOG_EVERY_SEC")) {
    char* end = nullptr;
    const long v = std::strtol(e, &end, 10);
    if (end != e && v >= 1) {
      return static_cast<int>(v);
    }
  }
  return 30;
}

std::string format_rate(double rate) {
  std::ostringstream os;
  if (rate >= 1e6) {
    os << std::fixed << std::setprecision(2) << (rate / 1e6) << "e6";
  } else if (rate >= 1e3) {
    os << std::fixed << std::setprecision(1) << (rate / 1e3) << "e3";
  } else {
    os << std::fixed << std::setprecision(2) << rate;
  }
  return os.str();
}

/// Wall-clock progress logger; safe to call from multiple threads (serialized logs).
class ProgressTicker {
 public:
  explicit ProgressTicker(std::string phase)
      : phase_(std::move(phase)),
        interval_sec_(progress_log_interval_sec()),
        start_(std::chrono::steady_clock::now()),
        last_log_(start_) {}

  void start_log(const std::string& detail = {}) {
    std::ostringstream os;
    os << phase_ << ": start";
    if (!detail.empty()) {
      os << " " << detail;
    }
    log_info(os.str());
  }

  /// processed units; total=nullopt → no percent/ETA.
  void tick(std::uint64_t processed, std::optional<std::uint64_t> total,
            const std::string& extra = {}) {
    maybe_log(processed, total, extra, false);
  }

  void done(std::uint64_t processed, std::optional<std::uint64_t> total,
            const std::string& extra = {}) {
    maybe_log(processed, total, extra, true);
  }

 private:
  void maybe_log(std::uint64_t processed, std::optional<std::uint64_t> total,
                 const std::string& extra, bool force) {
    const auto now = std::chrono::steady_clock::now();
    std::lock_guard<std::mutex> lock(mu_);
    const auto since_last =
        std::chrono::duration_cast<std::chrono::seconds>(now - last_log_).count();
    if (!force && since_last < interval_sec_) {
      return;
    }
    last_log_ = now;
    const double elapsed =
        std::chrono::duration<double>(now - start_).count();
    const double rate = elapsed > 0.0 ? static_cast<double>(processed) / elapsed : 0.0;

    std::ostringstream os;
    os << phase_ << ": ";
    if (total.has_value() && *total > 0) {
      const double pct = 100.0 * static_cast<double>(processed) / static_cast<double>(*total);
      os << processed << "/" << *total << " (" << std::fixed << std::setprecision(1) << pct
         << "%)";
    } else {
      os << "processed=" << processed;
    }
    os << " elapsed=" << static_cast<std::int64_t>(elapsed) << "s"
       << " rate=" << format_rate(rate) << "/s";
    if (total.has_value() && *total > processed && rate > 0.0) {
      const double eta = static_cast<double>(*total - processed) / rate;
      os << " eta=" << static_cast<std::int64_t>(eta) << "s";
    }
    if (!extra.empty()) {
      os << " " << extra;
    }
    if (force) {
      os << " done";
    }
    log_info(os.str());
  }

  std::string phase_;
  int interval_sec_;
  std::chrono::steady_clock::time_point start_;
  std::chrono::steady_clock::time_point last_log_;
  std::mutex mu_;
};

struct HeapItem {
  std::uint64_t code{0};
  std::size_t src{0};
};

struct HeapCmp {
  bool operator()(const HeapItem& a, const HeapItem& b) const {
    if (a.code != b.code) {
      return a.code > b.code;
    }
    return a.src > b.src;
  }
};

/// Buffered writer of sorted uint64 codes (one scatter bucket).
struct CodeBucketWriter {
  fs::path path;
  std::vector<std::uint64_t> buf;
  std::size_t batch{0};
  bool created{false};

  void init(fs::path p, std::size_t batch_rows) {
    path = std::move(p);
    batch = std::max<std::size_t>(1, batch_rows);
    buf.reserve(batch);
    created = false;
  }

  Error flush() {
    if (buf.empty()) {
      return Error::success();
    }
    std::ofstream out;
    if (!created) {
      out.open(path, std::ios::binary | std::ios::trunc);
      created = true;
    } else {
      out.open(path, std::ios::binary | std::ios::app);
    }
    if (!out) {
      return Error::io_error("failed to open scatter bucket: " + path.string() + " (" +
                             std::strerror(errno) + ")");
    }
    out.write(reinterpret_cast<const char*>(buf.data()),
              static_cast<std::streamsize>(buf.size() * sizeof(std::uint64_t)));
    if (!out) {
      return Error::io_error("failed writing scatter bucket: " + path.string());
    }
    out.close();
    buf.clear();
    return Error::success();
  }

  Error push(std::uint64_t code) {
    buf.push_back(code);
    if (buf.size() >= batch) {
      return flush();
    }
    return Error::success();
  }

  Error close() { return flush(); }
};

struct CodeCursor {
  std::ifstream in;
  std::uint64_t code{0};
  bool has{false};
  std::size_t acc{0};

  Error open(const fs::path& path, std::size_t accession) {
    acc = accession;
    has = false;
    std::error_code ec;
    if (!fs::exists(path, ec) || fs::file_size(path, ec) == 0) {
      return Error::success();
    }
    in.open(path, std::ios::binary);
    if (!in) {
      return Error::io_error("failed to open scatter bucket: " + path.string());
    }
    return advance();
  }

  Error advance() {
    has = false;
    in.read(reinterpret_cast<char*>(&code), sizeof(code));
    if (in.eof() && in.gcount() == 0) {
      return Error::success();
    }
    if (!in || static_cast<std::size_t>(in.gcount()) != sizeof(code)) {
      return Error::io_error("corrupt scatter bucket");
    }
    has = true;
    return Error::success();
  }
};

struct RowCursor {
  std::ifstream in;
  std::size_t n_stripes{0};
  std::size_t record_bytes{0};
  std::uint64_t code{0};
  std::vector<std::uint64_t> words;
  bool has{false};

  Error open(const fs::path& path, std::size_t stripes) {
    n_stripes = stripes;
    record_bytes = sizeof(std::uint64_t) * (1 + n_stripes);
    words.assign(n_stripes, 0);
    has = false;
    std::error_code ec;
    if (!fs::exists(path, ec) || fs::file_size(path, ec) == 0) {
      return Error::success();
    }
    in.open(path, std::ios::binary);
    if (!in) {
      return Error::io_error("failed to open merge spill: " + path.string());
    }
    return advance();
  }

  Error advance() {
    has = false;
    std::vector<char> rec(record_bytes);
    in.read(rec.data(), static_cast<std::streamsize>(record_bytes));
    if (in.eof() && in.gcount() == 0) {
      return Error::success();
    }
    if (!in || static_cast<std::size_t>(in.gcount()) != record_bytes) {
      return Error::io_error("corrupt merge spill record");
    }
    std::memcpy(&code, rec.data(), sizeof(code));
    std::memcpy(words.data(), rec.data() + sizeof(code), n_stripes * sizeof(std::uint64_t));
    has = true;
    return Error::success();
  }
};

Error write_row_file(std::ofstream& out, std::uint64_t code, const std::vector<std::uint64_t>& words) {
  out.write(reinterpret_cast<const char*>(&code), sizeof(code));
  out.write(reinterpret_cast<const char*>(words.data()),
            static_cast<std::streamsize>(words.size() * sizeof(std::uint64_t)));
  if (!out) {
    return Error::io_error("failed writing merge spill row");
  }
  return Error::success();
}

template <typename Fn>
Error merge_row_cursors(std::vector<RowCursor>& cursors, Fn&& on_row) {
  std::priority_queue<HeapItem, std::vector<HeapItem>, HeapCmp> heap;
  for (std::size_t i = 0; i < cursors.size(); ++i) {
    if (cursors[i].has) {
      heap.push(HeapItem{cursors[i].code, i});
    }
  }
  std::vector<std::uint64_t> words;
  if (!cursors.empty()) {
    words.assign(cursors[0].n_stripes, 0);
  }
  while (!heap.empty()) {
    const std::uint64_t code = heap.top().code;
    std::fill(words.begin(), words.end(), 0ULL);
    while (!heap.empty() && heap.top().code == code) {
      const HeapItem item = heap.top();
      heap.pop();
      auto& c = cursors[item.src];
      for (std::size_t s = 0; s < words.size(); ++s) {
        words[s] |= c.words[s];
      }
      if (auto err = c.advance(); !err.ok()) {
        return err;
      }
      if (c.has) {
        heap.push(HeapItem{c.code, item.src});
      }
    }
    if (auto err = on_row(code, words); !err.ok()) {
      return err;
    }
  }
  return Error::success();
}

Error merge_codes_to_row_spill(const std::vector<fs::path>& paths,
                               const std::vector<std::size_t>& acc_indices, std::size_t n_stripes,
                               const fs::path& out_path) {
  std::vector<CodeCursor> cursors(paths.size());
  for (std::size_t i = 0; i < paths.size(); ++i) {
    if (auto err = cursors[i].open(paths[i], acc_indices[i]); !err.ok()) {
      return err;
    }
  }
  std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to create code merge spill: " + out_path.string());
  }
  std::priority_queue<HeapItem, std::vector<HeapItem>, HeapCmp> heap;
  for (std::size_t i = 0; i < cursors.size(); ++i) {
    if (cursors[i].has) {
      heap.push(HeapItem{cursors[i].code, i});
    }
  }
  std::vector<std::uint64_t> words(n_stripes, 0);
  while (!heap.empty()) {
    const std::uint64_t code = heap.top().code;
    std::fill(words.begin(), words.end(), 0ULL);
    while (!heap.empty() && heap.top().code == code) {
      const HeapItem item = heap.top();
      heap.pop();
      const std::size_t acc = cursors[item.src].acc;
      words[acc / kAccessionsPerStripe] |= (1ULL << (acc % kAccessionsPerStripe));
      if (auto err = cursors[item.src].advance(); !err.ok()) {
        return err;
      }
      if (cursors[item.src].has) {
        heap.push(HeapItem{cursors[item.src].code, item.src});
      }
    }
    if (auto err = write_row_file(out, code, words); !err.ok()) {
      return err;
    }
  }
  return Error::success();
}

Error reduce_row_spills(std::vector<fs::path>& spills, std::size_t n_stripes, std::size_t fd_budget,
                        const fs::path& work, int& spill_seq) {
  const std::size_t group = std::max<std::size_t>(2, fd_budget);
  while (spills.size() > group) {
    std::vector<fs::path> next;
    for (std::size_t begin = 0; begin < spills.size(); begin += group) {
      const std::size_t end = std::min(spills.size(), begin + group);
      std::vector<RowCursor> cursors(end - begin);
      for (std::size_t i = begin; i < end; ++i) {
        if (auto err = cursors[i - begin].open(spills[i], n_stripes); !err.ok()) {
          return err;
        }
      }
      const fs::path out = work / ("reduce_" + std::to_string(spill_seq++) + ".rows");
      std::ofstream outf(out, std::ios::binary | std::ios::trunc);
      if (!outf) {
        return Error::io_error("failed to create reduce spill: " + out.string());
      }
      if (auto err = merge_row_cursors(cursors, [&](std::uint64_t code,
                                                   const std::vector<std::uint64_t>& words) {
            return write_row_file(outf, code, words);
          });
          !err.ok()) {
        return err;
      }
      for (std::size_t i = begin; i < end; ++i) {
        std::error_code ec;
        fs::remove(spills[i], ec);
      }
      next.push_back(out);
    }
    spills.swap(next);
  }
  return Error::success();
}

Error emit_pattern_row(std::ofstream& pat, std::ofstream& map,
                       std::unordered_map<PatternKey, std::uint32_t, PatternKeyHash>& dedup,
                       std::uint64_t code, const std::vector<std::uint64_t>& words,
                       std::size_t n_stripes) {
  PatternKey key;
  key.words = words;
  std::uint32_t pid = 0;
  const auto it = dedup.find(key);
  if (it == dedup.end()) {
    if (dedup.size() >= std::numeric_limits<std::uint32_t>::max()) {
      return Error::invalid_argument(
          "pattern count exceeds uint32 in one shard; raise --memory-gb or use more partitions");
    }
    pid = static_cast<std::uint32_t>(dedup.size());
    dedup.emplace(std::move(key), pid);
    pat.write(reinterpret_cast<const char*>(words.data()),
              static_cast<std::streamsize>(n_stripes * sizeof(std::uint64_t)));
    if (!pat) {
      return Error::io_error("failed writing pattern store shard");
    }
  } else {
    pid = it->second;
  }
  KmerMapEntry ent;
  ent.kmer_code = code;
  ent.pattern_id = pid;
  ent.pad = 0;
  map.write(reinterpret_cast<const char*>(&ent), sizeof(ent));
  if (!map) {
    return Error::io_error("failed writing k-mer map shard");
  }
  return Error::success();
}

/// Merge code buckets for one partition → pat/map (hierarchical if N > fd_budget).
Error merge_dedup_partition(std::size_t part, std::size_t n_acc, std::size_t n_stripes,
                            std::size_t fd_budget, const fs::path& scatter_root,
                            const fs::path& work, const fs::path& pat_path, const fs::path& map_path,
                            std::uint64_t& out_patterns, std::uint64_t& out_kmers) {
  out_patterns = 0;
  out_kmers = 0;

  std::vector<fs::path> bucket_paths;
  std::vector<std::size_t> acc_indices;
  bucket_paths.reserve(n_acc);
  acc_indices.reserve(n_acc);
  for (std::size_t i = 0; i < n_acc; ++i) {
    const fs::path p =
        scatter_root / ("a" + std::to_string(i)) / ("p" + std::to_string(part) + ".bin");
    std::error_code ec;
    if (fs::exists(p, ec) && fs::file_size(p, ec) > 0) {
      bucket_paths.push_back(p);
      acc_indices.push_back(i);
    }
  }
  if (bucket_paths.empty()) {
    std::ofstream(pat_path, std::ios::binary | std::ios::trunc);
    std::ofstream(map_path, std::ios::binary | std::ios::trunc);
    return Error::success();
  }

  std::ofstream pat(pat_path, std::ios::binary | std::ios::trunc);
  std::ofstream map(map_path, std::ios::binary | std::ios::trunc);
  if (!pat || !map) {
    return Error::io_error("failed to open shard outputs");
  }
  std::unordered_map<PatternKey, std::uint32_t, PatternKeyHash> dedup;
  dedup.reserve(1u << 16);

  auto on_row = [&](std::uint64_t code, const std::vector<std::uint64_t>& words) -> Error {
    ++out_kmers;
    return emit_pattern_row(pat, map, dedup, code, words, n_stripes);
  };

  if (bucket_paths.size() <= fd_budget) {
    std::vector<CodeCursor> cursors(bucket_paths.size());
    for (std::size_t i = 0; i < bucket_paths.size(); ++i) {
      if (auto err = cursors[i].open(bucket_paths[i], acc_indices[i]); !err.ok()) {
        return err;
      }
    }
    std::priority_queue<HeapItem, std::vector<HeapItem>, HeapCmp> heap;
    for (std::size_t i = 0; i < cursors.size(); ++i) {
      if (cursors[i].has) {
        heap.push(HeapItem{cursors[i].code, i});
      }
    }
    std::vector<std::uint64_t> words(n_stripes, 0);
    while (!heap.empty()) {
      const std::uint64_t code = heap.top().code;
      std::fill(words.begin(), words.end(), 0ULL);
      while (!heap.empty() && heap.top().code == code) {
        const HeapItem item = heap.top();
        heap.pop();
        const std::size_t acc = cursors[item.src].acc;
        words[acc / kAccessionsPerStripe] |= (1ULL << (acc % kAccessionsPerStripe));
        if (auto err = cursors[item.src].advance(); !err.ok()) {
          return err;
        }
        if (cursors[item.src].has) {
          heap.push(HeapItem{cursors[item.src].code, item.src});
        }
      }
      if (auto err = on_row(code, words); !err.ok()) {
        return err;
      }
    }
  } else {
    const fs::path part_work = work / ("p" + std::to_string(part) + "_hier");
    std::error_code ec;
    fs::create_directories(part_work, ec);
    const std::size_t group = std::max<std::size_t>(2, fd_budget);
    std::vector<fs::path> spills;
    int spill_seq = 0;
    for (std::size_t begin = 0; begin < bucket_paths.size(); begin += group) {
      const std::size_t end = std::min(bucket_paths.size(), begin + group);
      std::vector<fs::path> paths(bucket_paths.begin() + static_cast<std::ptrdiff_t>(begin),
                                  bucket_paths.begin() + static_cast<std::ptrdiff_t>(end));
      std::vector<std::size_t> idxs(acc_indices.begin() + static_cast<std::ptrdiff_t>(begin),
                                    acc_indices.begin() + static_cast<std::ptrdiff_t>(end));
      const fs::path spill = part_work / ("g" + std::to_string(spill_seq++) + ".rows");
      if (auto err = merge_codes_to_row_spill(paths, idxs, n_stripes, spill); !err.ok()) {
        return err;
      }
      spills.push_back(spill);
    }
    if (auto err = reduce_row_spills(spills, n_stripes, fd_budget, part_work, spill_seq);
        !err.ok()) {
      return err;
    }
    std::vector<RowCursor> cursors(spills.size());
    for (std::size_t i = 0; i < spills.size(); ++i) {
      if (auto err = cursors[i].open(spills[i], n_stripes); !err.ok()) {
        return err;
      }
    }
    if (auto err = merge_row_cursors(cursors, on_row); !err.ok()) {
      return err;
    }
    fs::remove_all(part_work, ec);
  }

  out_patterns = dedup.size();
  return Error::success();
}

struct MapCursor {
  std::ifstream in;
  KmerMapEntry ent{};
  bool has{false};
  std::uint32_t pid_offset{0};

  Error open(const fs::path& path, std::uint32_t offset) {
    pid_offset = offset;
    in.open(path, std::ios::binary);
    if (!in) {
      return Error::io_error("failed to open map shard: " + path.string());
    }
    return advance();
  }

  Error advance() {
    has = false;
    in.read(reinterpret_cast<char*>(&ent), sizeof(ent));
    if (in.eof() && in.gcount() == 0) {
      return Error::success();
    }
    if (!in) {
      return Error::io_error("failed reading map shard entry");
    }
    ent.pattern_id += pid_offset;
    has = true;
    return Error::success();
  }
};

struct MapHeapItem {
  std::uint64_t code{0};
  std::size_t shard{0};
};

struct MapHeapCmp {
  bool operator()(const MapHeapItem& a, const MapHeapItem& b) const {
    if (a.code != b.code) {
      return a.code > b.code;
    }
    return a.shard > b.shard;
  }
};

Error merge_map_files(const std::vector<fs::path>& paths, const std::vector<std::uint32_t>& offsets,
                      const fs::path& out_path) {
  std::vector<MapCursor> cursors(paths.size());
  std::priority_queue<MapHeapItem, std::vector<MapHeapItem>, MapHeapCmp> heap;
  for (std::size_t i = 0; i < paths.size(); ++i) {
    if (auto err = cursors[i].open(paths[i], offsets[i]); !err.ok()) {
      return err;
    }
    if (cursors[i].has) {
      heap.push(MapHeapItem{cursors[i].ent.kmer_code, i});
    }
  }
  std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to open map merge output: " + out_path.string());
  }
  while (!heap.empty()) {
    const auto item = heap.top();
    heap.pop();
    auto& c = cursors[item.shard];
    out.write(reinterpret_cast<const char*>(&c.ent), sizeof(KmerMapEntry));
    if (!out) {
      return Error::io_error("failed writing merged map");
    }
    if (auto err = c.advance(); !err.ok()) {
      return err;
    }
    if (c.has) {
      heap.push(MapHeapItem{c.ent.kmer_code, item.shard});
    }
  }
  return Error::success();
}

Error merge_maps_wave(std::vector<fs::path> paths, std::vector<std::uint32_t> offsets,
                      std::size_t fd_budget, const fs::path& work, std::ofstream& out,
                      ProgressTicker& progress, std::uint64_t total_kmers) {
  const std::size_t group = std::max<std::size_t>(2, fd_budget);
  int seq = 0;
  while (paths.size() > group) {
    std::vector<fs::path> next_paths;
    std::vector<std::uint32_t> next_offsets;
    for (std::size_t begin = 0; begin < paths.size(); begin += group) {
      const std::size_t end = std::min(paths.size(), begin + group);
      std::vector<fs::path> chunk(paths.begin() + static_cast<std::ptrdiff_t>(begin),
                                  paths.begin() + static_cast<std::ptrdiff_t>(end));
      std::vector<std::uint32_t> chunk_off(offsets.begin() + static_cast<std::ptrdiff_t>(begin),
                                           offsets.begin() + static_cast<std::ptrdiff_t>(end));
      const fs::path tmp = work / ("map_merge_" + std::to_string(seq++) + ".map");
      if (auto err = merge_map_files(chunk, chunk_off, tmp); !err.ok()) {
        return err;
      }
      for (std::size_t i = begin; i < end; ++i) {
        std::error_code ec;
        fs::remove(paths[i], ec);
      }
      next_paths.push_back(tmp);
      next_offsets.push_back(0);
    }
    paths.swap(next_paths);
    offsets.swap(next_offsets);
  }

  std::vector<MapCursor> cursors(paths.size());
  std::priority_queue<MapHeapItem, std::vector<MapHeapItem>, MapHeapCmp> heap;
  for (std::size_t i = 0; i < paths.size(); ++i) {
    if (auto err = cursors[i].open(paths[i], offsets[i]); !err.ok()) {
      return err;
    }
    if (cursors[i].has) {
      heap.push(MapHeapItem{cursors[i].ent.kmer_code, i});
    }
  }
  std::uint64_t written = 0;
  while (!heap.empty()) {
    const auto item = heap.top();
    heap.pop();
    auto& c = cursors[item.shard];
    out.write(reinterpret_cast<const char*>(&c.ent), sizeof(KmerMapEntry));
    if (!out) {
      return Error::io_error("failed writing k-mer map to matrix");
    }
    ++written;
    if ((written & 0xFFFFFull) == 0) {  // every ~1M rows, also gated by interval
      progress.tick(written, total_kmers);
    }
    if (auto err = c.advance(); !err.ok()) {
      return err;
    }
    if (c.has) {
      heap.push(MapHeapItem{c.ent.kmer_code, item.shard});
    }
  }
  progress.done(written, total_kmers);
  return Error::success();
}

Error assemble_v2(const BuildOptions& opts, const fs::path& work, std::size_t n_stripes,
                  std::size_t num_partitions, const std::vector<std::uint64_t>& pat_counts,
                  const std::vector<std::uint64_t>& kmer_counts, std::size_t fd_budget) {
  std::uint64_t total_patterns = 0;
  std::uint64_t total_kmers = 0;
  std::vector<std::uint32_t> offsets(num_partitions, 0);
  for (std::size_t i = 0; i < num_partitions; ++i) {
    if (total_patterns > std::numeric_limits<std::uint32_t>::max()) {
      return Error::invalid_argument("total pattern count exceeds uint32");
    }
    offsets[i] = static_cast<std::uint32_t>(total_patterns);
    total_patterns += pat_counts[i];
    total_kmers += kmer_counts[i];
  }
  if (total_patterns > std::numeric_limits<std::uint32_t>::max()) {
    return Error::invalid_argument("total pattern count exceeds uint32");
  }

  MatrixHeader header{};
  header.magic[0] = 'K';
  header.magic[1] = 'M';
  header.magic[2] = 'A';
  header.magic[3] = 'T';
  header.version = 2;
  header.kmer_size = static_cast<std::uint32_t>(opts.kmer_size);
  header.num_accessions = static_cast<std::uint32_t>(opts.accession_paths.size());
  header.num_stripes = static_cast<std::uint32_t>(n_stripes);
  header.num_rows = total_kmers;
  header.reserved = total_patterns;

  std::ofstream out(opts.output_path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to open matrix for writing: " + opts.output_path);
  }
  const std::size_t buf_sz = std::max<std::size_t>(1u << 16, runtime_config().io_buffer_bytes);
  std::vector<char> iobuf(buf_sz);
  out.rdbuf()->pubsetbuf(iobuf.data(), static_cast<std::streamsize>(iobuf.size()));

  out.write(reinterpret_cast<const char*>(&header), sizeof(header));
  if (!out) {
    return Error::io_error("failed writing matrix header");
  }

  ProgressTicker pat_progress("assemble_patterns");
  pat_progress.start_log("partitions=" + std::to_string(num_partitions));
  const std::size_t batch_rows = opts.batch_rows > 0 ? opts.batch_rows : 100000;
  std::vector<char> copy_buf(std::max<std::size_t>(1u << 20, batch_rows * n_stripes * 8));
  std::uint64_t pats_done = 0;
  for (std::size_t i = 0; i < num_partitions; ++i) {
    if (pat_counts[i] == 0) {
      ++pats_done;
      continue;
    }
    const fs::path pat = work / ("shard_" + std::to_string(i) + ".pat");
    std::ifstream in(pat, std::ios::binary);
    if (!in) {
      return Error::io_error("failed to open pattern shard: " + pat.string());
    }
    while (in) {
      in.read(copy_buf.data(), static_cast<std::streamsize>(copy_buf.size()));
      const auto n = in.gcount();
      if (n > 0) {
        out.write(copy_buf.data(), n);
        if (!out) {
          return Error::io_error("failed writing patterns to matrix");
        }
      }
    }
    ++pats_done;
    pat_progress.tick(pats_done, num_partitions);
  }
  pat_progress.done(pats_done, num_partitions);

  std::vector<fs::path> map_paths;
  std::vector<std::uint32_t> map_offsets;
  map_paths.reserve(num_partitions);
  for (std::size_t i = 0; i < num_partitions; ++i) {
    if (kmer_counts[i] == 0) {
      continue;
    }
    map_paths.push_back(work / ("shard_" + std::to_string(i) + ".map"));
    map_offsets.push_back(offsets[i]);
  }

  ProgressTicker map_progress("assemble");
  map_progress.start_log("map_rows_total=" + std::to_string(total_kmers));
  if (auto err = merge_maps_wave(std::move(map_paths), std::move(map_offsets), fd_budget, work, out,
                                 map_progress, total_kmers);
      !err.ok()) {
    return err;
  }

  return Error::success();
}

Error scatter_accession(std::size_t acc, const std::string& path, std::size_t num_partitions,
                        std::size_t batch, const fs::path& scatter_root, std::uint64_t& codes_out) {
  codes_out = 0;
  PresenceSetCursor cur;
  if (auto err = cur.open(path); !err.ok()) {
    return err;
  }
  const fs::path adir = scatter_root / ("a" + std::to_string(acc));
  std::error_code ec;
  fs::create_directories(adir, ec);
  if (ec) {
    return Error::io_error("failed to create scatter dir: " + adir.string());
  }

  std::vector<CodeBucketWriter> writers(num_partitions);
  for (std::size_t t = 0; t < num_partitions; ++t) {
    writers[t].init(adir / ("p" + std::to_string(t) + ".bin"), batch);
  }
  while (cur.has_value()) {
    const std::uint64_t code = cur.value();
    const std::size_t t = code_partition(code, num_partitions);
    if (auto err = writers[t].push(code); !err.ok()) {
      return err;
    }
    ++codes_out;
    if (auto err = cur.advance(); !err.ok()) {
      return err;
    }
  }
  for (std::size_t t = 0; t < num_partitions; ++t) {
    if (auto err = writers[t].close(); !err.ok()) {
      return err;
    }
    // Drop empty bucket files to cut inode count.
    if (!writers[t].created) {
      fs::remove(writers[t].path, ec);
    }
  }
  cur.close();
  return Error::success();
}

}  // namespace

Error build_matrix_from_presence_sets_streaming(const BuildOptions& opts) {
  if (opts.accession_paths.empty()) {
    return Error::invalid_argument("accession list is empty");
  }
  if (opts.kmer_size == 0 || opts.kmer_size > 32) {
    return Error::invalid_argument("k-mer size must be in 1..32");
  }
  if (opts.output_path.empty()) {
    return Error::invalid_argument("output path is required");
  }

  const std::size_t n_acc = opts.accession_paths.size();
  const std::size_t n_stripes = stripe_count_for_accessions(n_acc);
  const std::size_t threads =
      opts.num_threads > 0 ? opts.num_threads : effective_threads(runtime_config());
  const std::size_t memory_bytes =
      opts.memory_bytes > 0 ? opts.memory_bytes : default_memory_bytes();
  const std::size_t batch_rows_req = opts.batch_rows > 0 ? opts.batch_rows : 100000;
  const std::size_t fd_budget = merge_fd_budget();

  const fs::path work =
      fs::path(scratch_root(opts)) / ("kmat_build_" + std::to_string(::getpid()));
  std::error_code ec;
  fs::create_directories(work, ec);
  if (ec) {
    return Error::io_error("failed to create build work dir: " + work.string());
  }

  struct Cleaner {
    fs::path path;
    ~Cleaner() {
      std::error_code ignore;
      fs::remove_all(path, ignore);
    }
  } cleaner{work};

  std::uint64_t sum_kmers = 0;
  for (std::size_t i = 0; i < n_acc; ++i) {
    PresenceSetCursor cur;
    if (auto err = cur.open(opts.accession_paths[i]); !err.ok()) {
      return err;
    }
    if (cur.header().kmer_size != opts.kmer_size) {
      return Error::invalid_argument("presence set k-mer size mismatch: " +
                                     opts.accession_paths[i]);
    }
    sum_kmers += cur.header().num_kmers;
    cur.close();
  }

  const std::size_t mem_partitions =
      choose_partitions(n_stripes, memory_bytes, threads, batch_rows_req, sum_kmers);
  const std::size_t num_partitions = std::max(threads, mem_partitions);
  const std::size_t scatter_batch =
      scatter_batch_rows(batch_rows_req, num_partitions, threads, memory_bytes);

  {
    std::ostringstream msg;
    msg << "streaming build: accessions=" << n_acc << " stripes=" << n_stripes
        << " sum_kmers=" << sum_kmers << " partitions=" << num_partitions
        << " memory_gb=" << (memory_bytes / (1ull << 30)) << " scatter_batch=" << scatter_batch
        << " threads=" << threads << " merge_fd_budget=" << fd_budget
        << " log_every_sec=" << progress_log_interval_sec();
    log_info(msg.str());
  }

  const fs::path scatter_root = work / "scatter";
  fs::create_directories(scatter_root, ec);

  // --- Phase: scatter (parallel over accessions) ---
  ProgressTicker scatter_progress("scatter");
  scatter_progress.start_log("accessions=" + std::to_string(n_acc) +
                             " partitions=" + std::to_string(num_partitions) +
                             " threads=" + std::to_string(threads));
  std::atomic<std::uint64_t> accessions_done{0};
  std::atomic<std::uint64_t> codes_written{0};
  std::vector<Error> scatter_errors(n_acc, Error::success());

  parallel_for(0, n_acc, threads, [&](std::size_t i) {
    std::uint64_t codes = 0;
    if (auto err =
            scatter_accession(i, opts.accession_paths[i], num_partitions, scatter_batch,
                              scatter_root, codes);
        !err.ok()) {
      scatter_errors[i] = err;
      return;
    }
    codes_written.fetch_add(codes, std::memory_order_relaxed);
    const std::uint64_t done = accessions_done.fetch_add(1, std::memory_order_relaxed) + 1;
    scatter_progress.tick(done, n_acc, "codes=" + std::to_string(codes_written.load()));
  });

  for (const Error& err : scatter_errors) {
    if (!err.ok()) {
      return err;
    }
  }
  scatter_progress.done(n_acc, n_acc, "codes=" + std::to_string(codes_written.load()));

  // --- Phase: parallel merge + pattern dedup per code partition ---
  const std::size_t merge_workers =
      std::min(threads, std::max<std::size_t>(1, fd_budget / 4));
  ProgressTicker merge_progress("merge_dedup");
  merge_progress.start_log("partitions=" + std::to_string(num_partitions) +
                           " workers=" + std::to_string(merge_workers));

  std::vector<std::uint64_t> pat_counts(num_partitions, 0);
  std::vector<std::uint64_t> kmer_counts(num_partitions, 0);
  std::vector<Error> merge_errors(num_partitions, Error::success());
  std::atomic<std::uint64_t> parts_done{0};

  parallel_for(0, num_partitions, merge_workers, [&](std::size_t p) {
    const fs::path pat = work / ("shard_" + std::to_string(p) + ".pat");
    const fs::path map = work / ("shard_" + std::to_string(p) + ".map");
    if (auto err = merge_dedup_partition(p, n_acc, n_stripes, fd_budget, scatter_root, work, pat,
                                         map, pat_counts[p], kmer_counts[p]);
        !err.ok()) {
      merge_errors[p] = err;
    }
    const std::uint64_t done = parts_done.fetch_add(1, std::memory_order_relaxed) + 1;
    merge_progress.tick(done, num_partitions);
  });

  for (const Error& err : merge_errors) {
    if (!err.ok()) {
      return err;
    }
  }
  const std::uint64_t unique_emitted =
      std::accumulate(kmer_counts.begin(), kmer_counts.end(), 0ull);
  merge_progress.done(num_partitions, num_partitions,
                      "unique_kmers=" + std::to_string(unique_emitted));

  // Scatter trees no longer needed.
  fs::remove_all(scatter_root, ec);

  log_info("assembling v2 matrix: patterns=" +
           std::to_string(std::accumulate(pat_counts.begin(), pat_counts.end(), 0ull)) +
           " kmers=" + std::to_string(unique_emitted));

  return assemble_v2(opts, work, n_stripes, num_partitions, pat_counts, kmer_counts, fd_budget);
}

}  // namespace kmat
