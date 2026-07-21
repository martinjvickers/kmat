#include "kmat/universe.hpp"

#include "kmat/log.hpp"
#include "kmat/presence.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <filesystem>
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

Error merge_sorted_code_streams(std::vector<PresenceSetCursor>& cursors, std::ofstream& out,
                                std::uint64_t& unique_out) {
  unique_out = 0;
  std::priority_queue<HeapItem, std::vector<HeapItem>, HeapCmp> heap;
  for (std::size_t i = 0; i < cursors.size(); ++i) {
    if (cursors[i].has_value()) {
      heap.push(HeapItem{cursors[i].value(), i});
    }
  }
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
    out.write(reinterpret_cast<const char*>(&code), sizeof(code));
    if (!out) {
      return Error::io_error("failed writing universe codes");
    }
    ++unique_out;
  }
  return Error::success();
}

Error merge_kset_group(const std::vector<std::string>& paths, std::size_t kmer_size,
                       const fs::path& out_path) {
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
  header.num_kmers = 0;  // rewrite after count

  std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to create universe: " + out_path.string());
  }
  out.write(reinterpret_cast<const char*>(&header), sizeof(header));
  if (!out) {
    return Error::io_error("failed writing universe header");
  }

  std::uint64_t unique = 0;
  if (auto err = merge_sorted_code_streams(cursors, out, unique); !err.ok()) {
    return err;
  }
  for (auto& c : cursors) {
    c.close();
  }

  header.num_kmers = unique;
  out.seekp(0);
  out.write(reinterpret_cast<const char*>(&header), sizeof(header));
  if (!out) {
    return Error::io_error("failed rewriting universe header count");
  }
  return Error::success();
}

Error merge_kuniv_group(const std::vector<fs::path>& paths, std::size_t kmer_size,
                        const fs::path& out_path) {
  struct UCursor {
    UniverseCursor cur;
    bool open_ok{false};
  };
  std::vector<UCursor> cursors(paths.size());
  for (std::size_t i = 0; i < paths.size(); ++i) {
    if (auto err = cursors[i].cur.open(paths[i].string()); !err.ok()) {
      return err;
    }
    if (cursors[i].cur.header().kmer_size != kmer_size) {
      return Error::invalid_argument("k-mer size mismatch in universe spill");
    }
    cursors[i].open_ok = true;
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
    out.write(reinterpret_cast<const char*>(&code), sizeof(code));
    if (!out) {
      return Error::io_error("failed writing reduced universe");
    }
    ++unique;
  }
  for (auto& c : cursors) {
    c.cur.close();
  }
  header.num_kmers = unique;
  out.seekp(0);
  out.write(reinterpret_cast<const char*>(&header), sizeof(header));
  return Error::success();
}

}  // namespace

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

  log_info("build-master: accessions=" + std::to_string(kset_paths.size()) +
           " group_size=" + std::to_string(G));

  const auto t0 = std::chrono::steady_clock::now();
  std::vector<fs::path> level;
  int seq = 0;
  for (std::size_t begin = 0; begin < kset_paths.size(); begin += G) {
    const std::size_t end = std::min(kset_paths.size(), begin + G);
    std::vector<std::string> group(kset_paths.begin() + static_cast<std::ptrdiff_t>(begin),
                                   kset_paths.begin() + static_cast<std::ptrdiff_t>(end));
    const fs::path out = work / ("m" + std::to_string(seq++) + ".kuniv");
    if (auto err = merge_kset_group(group, kmer_size, out); !err.ok()) {
      return err;
    }
    level.push_back(out);
    const auto elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0)
            .count();
    log_info("build-master: groups=" + std::to_string(level.size()) + "/" +
             std::to_string((kset_paths.size() + G - 1) / G) + " elapsed=" +
             std::to_string(elapsed) + "s");
  }

  while (level.size() > 1) {
    std::vector<fs::path> next;
    for (std::size_t begin = 0; begin < level.size(); begin += G) {
      const std::size_t end = std::min(level.size(), begin + G);
      std::vector<fs::path> group(level.begin() + static_cast<std::ptrdiff_t>(begin),
                                  level.begin() + static_cast<std::ptrdiff_t>(end));
      const fs::path out = work / ("r" + std::to_string(seq++) + ".kuniv");
      if (auto err = merge_kuniv_group(group, kmer_size, out); !err.ok()) {
        return err;
      }
      for (const auto& p : group) {
        fs::remove(p, ec);
      }
      next.push_back(out);
    }
    level.swap(next);
    log_info("build-master: reduce level size=" + std::to_string(level.size()));
  }

  fs::copy_file(level.front(), output_path, fs::copy_options::overwrite_existing, ec);
  if (ec) {
    // Fallback: rename if same filesystem
    fs::rename(level.front(), output_path, ec);
    if (ec) {
      return Error::io_error("failed to write universe output: " + output_path);
    }
  }

  UniverseHeader hdr{};
  if (auto err = read_universe_header(output_path, hdr); !err.ok()) {
    return err;
  }
  log_info("build-master: done unique_kmers=" + std::to_string(hdr.num_kmers) + " path=" +
           output_path);
  return Error::success();
}

}  // namespace kmat
