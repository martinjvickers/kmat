#include "kmat/universe.hpp"

#include "kmat/log.hpp"
#include "kmat/presence.hpp"
#include "kmat/runtime.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <queue>
#include <sstream>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

namespace kmat {

namespace {

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

std::string scratch_root(const std::string& tmpdir) {
  if (!tmpdir.empty()) {
    return tmpdir;
  }
  if (const char* t = std::getenv("TMPDIR")) {
    if (*t) {
      return t;
    }
  }
  return "/tmp";
}

class MergeProgress {
 public:
  explicit MergeProgress(std::string label) : label_(std::move(label)), t0_(Clock::now()) {}

  void tick(std::uint64_t unique_emitted) {
    const auto now = Clock::now();
    if (now - last_ < std::chrono::seconds(15) && unique_emitted - last_unique_ < 5'000'000) {
      return;
    }
    last_ = now;
    last_unique_ = unique_emitted;
    const auto sec = std::chrono::duration_cast<std::chrono::seconds>(now - t0_).count();
    const double rate = sec > 0 ? static_cast<double>(unique_emitted) / static_cast<double>(sec) : 0.0;
    std::ostringstream oss;
    oss << label_ << ": unique=" << unique_emitted << " elapsed=" << sec << "s rate="
        << static_cast<std::uint64_t>(rate) << "/s";
    log_info(oss.str());
  }

  void done(std::uint64_t unique_emitted) {
    const auto sec =
        std::chrono::duration_cast<std::chrono::seconds>(Clock::now() - t0_).count();
    log_info(label_ + ": done unique=" + std::to_string(unique_emitted) + " elapsed=" +
             std::to_string(sec) + "s");
  }

 private:
  using Clock = std::chrono::steady_clock;
  std::string label_;
  Clock::time_point t0_;
  Clock::time_point last_{t0_};
  std::uint64_t last_unique_{0};
};

Error flush_code_buf(std::ofstream& out, std::vector<std::uint64_t>& buf) {
  if (buf.empty()) {
    return Error::success();
  }
  out.write(reinterpret_cast<const char*>(buf.data()),
            static_cast<std::streamsize>(buf.size() * sizeof(std::uint64_t)));
  if (!out) {
    return Error::io_error("failed writing universe codes");
  }
  buf.clear();
  return Error::success();
}

/// Patch num_kmers after the payload is fully closed. In-place seekp on a large
/// custom-buffered ofstream can leave num_kmers=0 on disk while the file body is huge.
Error rewrite_universe_num_kmers(const fs::path& path, std::uint64_t num_kmers) {
  std::fstream f(path, std::ios::binary | std::ios::in | std::ios::out);
  if (!f) {
    return Error::io_error("failed to reopen universe for header rewrite: " + path.string());
  }
  UniverseHeader header{};
  f.read(reinterpret_cast<char*>(&header), sizeof(header));
  if (!f) {
    return Error::io_error("failed reading universe header for rewrite: " + path.string());
  }
  if (header.magic[0] != 'K' || header.magic[1] != 'U' || header.magic[2] != 'N' ||
      header.magic[3] != 'I') {
    return Error::invalid_argument("invalid universe magic during header rewrite");
  }
  header.num_kmers = num_kmers;
  f.seekp(0);
  f.write(reinterpret_cast<const char*>(&header), sizeof(header));
  f.flush();
  if (!f) {
    return Error::io_error("failed rewriting universe num_kmers: " + path.string());
  }
  return Error::success();
}

Error merge_sorted_code_streams(std::vector<PresenceSetCursor>& cursors, std::ofstream& out,
                                std::uint64_t& unique_out, MergeProgress* progress) {
  unique_out = 0;
  std::priority_queue<HeapItem, std::vector<HeapItem>, HeapCmp> heap;
  for (std::size_t i = 0; i < cursors.size(); ++i) {
    if (cursors[i].has_value()) {
      heap.push(HeapItem{cursors[i].value(), i});
    }
  }
  std::vector<std::uint64_t> buf;
  buf.reserve(1u << 20);
  while (!heap.empty()) {
    const std::uint64_t code = heap.top().code;
    while (!heap.empty() && heap.top().code == code) {
      const HeapItem item = heap.top();
      heap.pop();
      if (auto err = cursors[item.src].advance(); !err.ok()) {
        return err;
      }
      if (cursors[item.src].has_value()) {
        heap.push(HeapItem{cursors[item.src].value(), item.src});
      }
    }
    buf.push_back(code);
    ++unique_out;
    if (buf.size() >= (1u << 20)) {
      if (auto err = flush_code_buf(out, buf); !err.ok()) {
        return err;
      }
      if (progress) {
        progress->tick(unique_out);
      }
    }
  }
  if (auto err = flush_code_buf(out, buf); !err.ok()) {
    return err;
  }
  return Error::success();
}

Error merge_kset_group(const std::vector<std::string>& paths, std::size_t kmer_size,
                       const fs::path& out_path, const std::string& progress_label) {
  if (paths.empty()) {
    return Error::invalid_argument("empty kset group");
  }
  std::vector<PresenceSetCursor> cursors(paths.size());
  for (std::size_t i = 0; i < paths.size(); ++i) {
    if (auto err = cursors[i].open(paths[i]); !err.ok()) {
      return err;
    }
    if (cursors[i].header().kmer_size != kmer_size) {
      return Error::invalid_argument("k-mer size mismatch in presence set: " + paths[i]);
    }
  }

  UniverseHeader header{};
  header.magic[0] = 'K';
  header.magic[1] = 'U';
  header.magic[2] = 'N';
  header.magic[3] = 'I';
  header.version = 1;
  header.kmer_size = static_cast<std::uint32_t>(kmer_size);
  header.num_kmers = 0;  // patched after close

  std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to create universe: " + out_path.string());
  }

  out.write(reinterpret_cast<const char*>(&header), sizeof(header));
  if (!out) {
    return Error::io_error("failed writing universe header");
  }

  MergeProgress progress(progress_label);
  std::uint64_t unique = 0;
  if (auto err = merge_sorted_code_streams(cursors, out, unique, &progress); !err.ok()) {
    return err;
  }
  for (auto& c : cursors) {
    c.close();
  }
  out.flush();
  out.close();
  if (auto err = rewrite_universe_num_kmers(out_path, unique); !err.ok()) {
    return err;
  }
  progress.done(unique);
  return Error::success();
}

Error merge_kuniv_group(const std::vector<fs::path>& paths, std::size_t kmer_size,
                        const fs::path& out_path, const std::string& progress_label) {
  struct UCursor {
    UniverseCursor cur;
  };
  std::vector<UCursor> cursors(paths.size());
  for (std::size_t i = 0; i < paths.size(); ++i) {
    if (auto err = cursors[i].cur.open(paths[i].string()); !err.ok()) {
      return err;
    }
    if (cursors[i].cur.header().kmer_size != kmer_size) {
      return Error::invalid_argument("k-mer size mismatch in universe spill");
    }
  }

  UniverseHeader header{};
  header.magic[0] = 'K';
  header.magic[1] = 'U';
  header.magic[2] = 'N';
  header.magic[3] = 'I';
  header.version = 1;
  header.kmer_size = static_cast<std::uint32_t>(kmer_size);

  std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to create universe reduce: " + out_path.string());
  }
  out.write(reinterpret_cast<const char*>(&header), sizeof(header));

  std::priority_queue<HeapItem, std::vector<HeapItem>, HeapCmp> heap;
  for (std::size_t i = 0; i < cursors.size(); ++i) {
    if (cursors[i].cur.has_value()) {
      heap.push(HeapItem{cursors[i].cur.value(), i});
    }
  }
  MergeProgress progress(progress_label);
  std::vector<std::uint64_t> buf;
  buf.reserve(1u << 20);
  std::uint64_t unique = 0;
  while (!heap.empty()) {
    const std::uint64_t code = heap.top().code;
    while (!heap.empty() && heap.top().code == code) {
      const HeapItem item = heap.top();
      heap.pop();
      if (auto err = cursors[item.src].cur.advance(); !err.ok()) {
        return err;
      }
      if (cursors[item.src].cur.has_value()) {
        heap.push(HeapItem{cursors[item.src].cur.value(), item.src});
      }
    }
    buf.push_back(code);
    ++unique;
    if (buf.size() >= (1u << 20)) {
      if (auto err = flush_code_buf(out, buf); !err.ok()) {
        return err;
      }
      progress.tick(unique);
    }
  }
  if (auto err = flush_code_buf(out, buf); !err.ok()) {
    return err;
  }
  for (auto& c : cursors) {
    c.cur.close();
  }
  out.flush();
  out.close();
  if (auto err = rewrite_universe_num_kmers(out_path, unique); !err.ok()) {
    return err;
  }
  progress.done(unique);
  return Error::success();
}

}  // namespace

Error payload_count_from_file(const std::string& path, std::uint64_t& count_out) {
  std::ifstream in(path, std::ios::binary | std::ios::ate);
  if (!in) {
    return Error::io_error("failed to stat universe: " + path);
  }
  const auto sz = static_cast<std::uint64_t>(in.tellg());
  if (sz < sizeof(UniverseHeader)) {
    return Error::invalid_argument("universe file too small: " + path);
  }
  const std::uint64_t payload = sz - sizeof(UniverseHeader);
  if (payload % sizeof(std::uint64_t) != 0) {
    return Error::invalid_argument("universe payload size not a multiple of 8: " + path);
  }
  count_out = payload / sizeof(std::uint64_t);
  return Error::success();
}

Error write_universe(const std::string& path, const UniverseSet& set) {
  if (set.header.num_kmers != set.kmers.size()) {
    return Error::invalid_argument("universe header count mismatch");
  }
  for (std::size_t i = 1; i < set.kmers.size(); ++i) {
    if (set.kmers[i] <= set.kmers[i - 1]) {
      return Error::invalid_argument("universe k-mers must be strictly sorted");
    }
  }
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to open universe for writing: " + path);
  }
  UniverseHeader header = set.header;
  header.magic[0] = 'K';
  header.magic[1] = 'U';
  header.magic[2] = 'N';
  header.magic[3] = 'I';
  header.version = 1;
  header.num_kmers = set.kmers.size();
  out.write(reinterpret_cast<const char*>(&header), sizeof(header));
  if (!set.kmers.empty()) {
    out.write(reinterpret_cast<const char*>(set.kmers.data()),
              static_cast<std::streamsize>(set.kmers.size() * sizeof(std::uint64_t)));
  }
  if (!out) {
    return Error::io_error("failed writing universe: " + path);
  }
  return Error::success();
}

Error read_universe_header(const std::string& path, UniverseHeader& header) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return Error::io_error("failed to open universe: " + path);
  }
  in.read(reinterpret_cast<char*>(&header), sizeof(header));
  if (!in) {
    return Error::io_error("failed reading universe header");
  }
  if (header.magic[0] != 'K' || header.magic[1] != 'U' || header.magic[2] != 'N' ||
      header.magic[3] != 'I') {
    return Error::invalid_argument("invalid universe magic (expected KUNI)");
  }
  if (header.version != 1) {
    return Error::invalid_argument("unsupported universe version");
  }
  std::uint64_t payload = 0;
  if (auto err = payload_count_from_file(path, payload); !err.ok()) {
    return err;
  }
  if (header.num_kmers == 0 && payload > 0) {
    log_warn("universe header num_kmers=0 but file has " + std::to_string(payload) +
             " codes; recovering from file size (" + path + ")");
    header.num_kmers = payload;
  } else if (header.num_kmers != payload) {
    return Error::invalid_argument("universe header count (" + std::to_string(header.num_kmers) +
                                   ") != file payload (" + std::to_string(payload) + "): " + path);
  }
  return Error::success();
}

Error UniverseCursor::open(const std::string& path) {
  close();
  in_.open(path, std::ios::binary);
  if (!in_) {
    return Error::io_error("failed to open universe: " + path);
  }
  in_.read(reinterpret_cast<char*>(&header_), sizeof(header_));
  if (!in_) {
    return Error::io_error("failed reading universe header");
  }
  if (header_.magic[0] != 'K' || header_.magic[1] != 'U' || header_.magic[2] != 'N' ||
      header_.magic[3] != 'I') {
    return Error::invalid_argument("invalid universe magic");
  }
  if (header_.version != 1) {
    return Error::invalid_argument("unsupported universe version");
  }
  std::uint64_t payload = 0;
  if (auto err = payload_count_from_file(path, payload); !err.ok()) {
    return err;
  }
  if (header_.num_kmers == 0 && payload > 0) {
    log_warn("universe header num_kmers=0 but file has " + std::to_string(payload) +
             " codes; recovering from file size (" + path + ")");
    header_.num_kmers = payload;
  } else if (header_.num_kmers != payload) {
    return Error::invalid_argument("universe header count (" + std::to_string(header_.num_kmers) +
                                   ") != file payload (" + std::to_string(payload) + "): " + path);
  }
  remaining_ = header_.num_kmers;
  return advance();
}

void UniverseCursor::close() {
  if (in_.is_open()) {
    in_.close();
  }
  has_value_ = false;
  remaining_ = 0;
}

Error UniverseCursor::advance() {
  has_value_ = false;
  if (remaining_ == 0) {
    return Error::success();
  }
  in_.read(reinterpret_cast<char*>(&value_), sizeof(value_));
  if (!in_) {
    return Error::io_error("failed reading universe code");
  }
  --remaining_;
  has_value_ = true;
  return Error::success();
}

bool path_looks_universe(const std::string& path) {
  if (path.size() >= 6) {
    const auto ext = path.substr(path.size() - 6);
    if (ext == ".kuniv" || ext == ".KUNIV") {
      return true;
    }
  }
  return false;
}

Error build_universe_from_presence_sets(const std::vector<std::string>& kset_paths,
                                        std::size_t kmer_size, const std::string& output_path,
                                        const std::string& tmpdir, std::size_t group_size) {
  if (kset_paths.empty()) {
    return Error::invalid_argument("accession list is empty");
  }
  if (kmer_size == 0 || kmer_size > 32) {
    return Error::invalid_argument("k-mer size must be in 1..32");
  }
  const std::size_t G = std::max<std::size_t>(2, group_size);
  const std::size_t threads = effective_threads(runtime_config());

  const fs::path work =
      fs::path(scratch_root(tmpdir)) / ("kmat_kuniv_" + std::to_string(::getpid()));
  std::error_code ec;
  fs::create_directories(work, ec);
  if (ec) {
    return Error::io_error("failed to create universe work dir: " + work.string());
  }
  struct Cleaner {
    fs::path path;
    ~Cleaner() {
      std::error_code ignore;
      fs::remove_all(path, ignore);
    }
  } cleaner{work};

  const std::size_t num_groups = (kset_paths.size() + G - 1) / G;
  log_info("build-master: accessions=" + std::to_string(kset_paths.size()) +
           " groups=" + std::to_string(num_groups) + " group_size=" + std::to_string(G) +
           " threads=" + std::to_string(threads));

  const auto t0 = std::chrono::steady_clock::now();
  std::vector<fs::path> level(num_groups);
  std::vector<Error> group_errs(num_groups, Error::success());
  std::atomic<std::size_t> groups_done{0};
  std::mutex log_mu;

  // Independent group merges — this is where most wall-clock time is for large panels.
  parallel_for(0, num_groups, threads, [&](std::size_t gi) {
    const std::size_t begin = gi * G;
    const std::size_t end = std::min(kset_paths.size(), begin + G);
    std::vector<std::string> group(kset_paths.begin() + static_cast<std::ptrdiff_t>(begin),
                                   kset_paths.begin() + static_cast<std::ptrdiff_t>(end));
    const fs::path out = work / ("m" + std::to_string(gi) + ".kuniv");
    const std::string label = "build-master group " + std::to_string(gi + 1) + "/" +
                              std::to_string(num_groups) + " (n=" + std::to_string(group.size()) +
                              ")";
    group_errs[gi] = merge_kset_group(group, kmer_size, out, label);
    if (!group_errs[gi].ok()) {
      return;
    }
    level[gi] = out;
    const std::size_t done = groups_done.fetch_add(1) + 1;
    const auto elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0)
            .count();
    std::lock_guard<std::mutex> lock(log_mu);
    log_info("build-master: groups_done=" + std::to_string(done) + "/" +
             std::to_string(num_groups) + " elapsed=" + std::to_string(elapsed) + "s");
  });

  for (const auto& err : group_errs) {
    if (!err.ok()) {
      return err;
    }
  }

  int seq = static_cast<int>(num_groups);
  while (level.size() > 1) {
    const std::size_t n_in = level.size();
    const std::size_t n_out = (n_in + G - 1) / G;
    log_info("build-master: reduce level in=" + std::to_string(n_in) + " out=" +
             std::to_string(n_out) + " threads=" + std::to_string(std::min(threads, n_out)));

    std::vector<fs::path> next(n_out);
    std::vector<Error> reduce_errs(n_out, Error::success());
    const int seq_base = seq;
    seq += static_cast<int>(n_out);

    parallel_for(0, n_out, std::min(threads, n_out), [&](std::size_t ri) {
      const std::size_t begin = ri * G;
      const std::size_t end = std::min(n_in, begin + G);
      std::vector<fs::path> group(level.begin() + static_cast<std::ptrdiff_t>(begin),
                                  level.begin() + static_cast<std::ptrdiff_t>(end));
      const fs::path out = work / ("r" + std::to_string(seq_base + static_cast<int>(ri)) + ".kuniv");
      const std::string label = "build-master reduce " + std::to_string(ri + 1) + "/" +
                                std::to_string(n_out);
      reduce_errs[ri] = merge_kuniv_group(group, kmer_size, out, label);
      if (!reduce_errs[ri].ok()) {
        return;
      }
      next[ri] = out;
      for (const auto& p : group) {
        std::error_code ignore;
        fs::remove(p, ignore);
      }
    });

    for (const auto& err : reduce_errs) {
      if (!err.ok()) {
        return err;
      }
    }
    level.swap(next);
  }

  fs::copy_file(level.front(), output_path, fs::copy_options::overwrite_existing, ec);
  if (ec) {
    fs::rename(level.front(), output_path, ec);
    if (ec) {
      return Error::io_error("failed to write universe output: " + output_path);
    }
  }

  UniverseHeader hdr{};
  if (auto err = read_universe_header(output_path, hdr); !err.ok()) {
    return err;
  }
  if (hdr.num_kmers == 0) {
    return Error::invalid_argument(
        "build-master produced an empty universe (0 unique k-mers); check that .kset inputs "
        "are non-empty");
  }
  const auto elapsed =
      std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0)
          .count();
  log_info("build-master: done unique_kmers=" + std::to_string(hdr.num_kmers) + " elapsed=" +
           std::to_string(elapsed) + "s path=" + output_path);
  return Error::success();
}

}  // namespace kmat
