#include "kmat/build_stream.hpp"

#include "kmat/log.hpp"
#include "kmat/matrix_layout.hpp"
#include "kmat/presence.hpp"
#include "kmat/runtime.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
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

std::size_t hash_words(const std::vector<std::uint64_t>& words) {
  PatternKey key;
  key.words = words;
  return PatternKeyHash{}(key);
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
    return 64ull << 30;  // 64 GiB
  }
  return 8ull << 30;  // 8 GiB
}

/// Raise soft RLIMIT_NOFILE to the hard limit (common on HPC where soft=1024).
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

/// How many files we may hold open for a multiway merge (exclusive of stdin/out/err).
/// Override with KMAT_BUILD_MAX_OPEN for tests / forced hierarchical merge.
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

/// Partition count from memory budget. ΣK_i is only an *upper bound* on unique U —
/// using it makes P larger (safer for RAM), never opens that many FDs (shard writers
/// open-append-close; map assemble merges in waves if P is large).
std::size_t choose_partitions(std::size_t n_stripes, std::size_t memory_bytes,
                              std::size_t num_threads, std::size_t batch_rows,
                              std::uint64_t sum_kmers) {
  const std::size_t threads = std::max<std::size_t>(1, num_threads);
  const std::size_t record_bytes = sizeof(std::uint64_t) * (1 + n_stripes);
  const std::size_t batch = std::max<std::size_t>(1, batch_rows);

  // Target: expected rows/shard ≈ U/P fit in ~memory/threads (hash table ~2× row bytes).
  const std::size_t budget_per_thread =
      std::max<std::size_t>(1ull << 20, memory_bytes / threads);
  const std::uint64_t u_upper = std::max<std::uint64_t>(1, sum_kmers);
  const std::uint64_t bytes_per_row_est = record_bytes * 2;
  const std::uint64_t rows_per_shard =
      std::max<std::uint64_t>(1, budget_per_thread / bytes_per_row_est);
  std::uint64_t p64 = (u_upper + rows_per_shard - 1) / rows_per_shard;
  p64 = std::max<std::uint64_t>(threads, p64);

  // Cap by aggregate writer buffer budget (P × batch × record), not by FD count.
  const std::size_t buf_budget = std::max<std::size_t>(record_bytes * batch, memory_bytes / 4);
  const std::uint64_t max_p_by_buf =
      std::max<std::uint64_t>(1, buf_budget / (record_bytes * batch));
  p64 = std::min(p64, max_p_by_buf);

  // Absolute ceiling: enough for huge panels; still far below inode/path pain.
  constexpr std::uint64_t kAbsMaxPartitions = 1ull << 16;  // 65536
  p64 = std::min(p64, kAbsMaxPartitions);
  p64 = std::max<std::uint64_t>(1, p64);
  return static_cast<std::size_t>(p64);
}

/// Shrink per-shard flush quantum so P writer buffers stay inside ~memory/4.
std::size_t effective_batch_rows(std::size_t batch_rows, std::size_t num_partitions,
                                 std::size_t n_stripes, std::size_t memory_bytes) {
  const std::size_t record_bytes = sizeof(std::uint64_t) * (1 + n_stripes);
  const std::size_t batch = std::max<std::size_t>(1, batch_rows);
  const std::size_t p = std::max<std::size_t>(1, num_partitions);
  const std::size_t buf_budget = std::max<std::size_t>(record_bytes, memory_bytes / 4);
  const std::size_t max_batch = std::max<std::size_t>(1, buf_budget / (p * record_bytes));
  return std::min(batch, max_batch);
}

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

struct ShardWriter {
  fs::path path;
  std::vector<char> buf;
  std::size_t buffered_rows{0};
  std::size_t batch_rows{0};
  std::size_t record_bytes{0};
  std::uint64_t total_rows{0};
  bool created{false};

  Error init(const fs::path& shard_path, std::size_t n_stripes, std::size_t batch) {
    path = shard_path;
    record_bytes = sizeof(std::uint64_t) * (1 + n_stripes);
    batch_rows = std::max<std::size_t>(1, batch);
    buf.reserve(batch_rows * record_bytes);
    created = false;
    return Error::success();
  }

  /// Open-append-write-close: partition count is not limited by FDs.
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
      return Error::io_error("failed to open shard for writing: " + path.string() + " (" +
                             std::strerror(errno) + ")");
    }
    out.write(buf.data(), static_cast<std::streamsize>(buf.size()));
    if (!out) {
      return Error::io_error("failed writing shard rows: " + path.string());
    }
    out.close();
    buf.clear();
    buffered_rows = 0;
    return Error::success();
  }

  Error write_row(std::uint64_t code, const std::vector<std::uint64_t>& words) {
    const std::size_t off = buf.size();
    buf.resize(off + record_bytes);
    std::memcpy(buf.data() + off, &code, sizeof(code));
    std::memcpy(buf.data() + off + sizeof(code), words.data(),
                words.size() * sizeof(std::uint64_t));
    ++buffered_rows;
    ++total_rows;
    if (buffered_rows >= batch_rows) {
      return flush();
    }
    return Error::success();
  }

  Error close() { return flush(); }
};

/// Sorted stream of (code, stripe-words) rows — used for hierarchical merge spills.
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

/// Multiway-merge sorted row streams; OR stripe words for equal codes; invoke on_row.
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

/// Merge a subset of .kset files into one sorted (code, words) spill.
Error merge_ksets_to_spill(const std::vector<std::string>& paths,
                           const std::vector<std::size_t>& acc_indices, std::size_t n_stripes,
                           const fs::path& out_path) {
  std::vector<PresenceSetCursor> cursors(paths.size());
  for (std::size_t i = 0; i < paths.size(); ++i) {
    if (auto err = cursors[i].open(paths[i]); !err.ok()) {
      return err;
    }
  }

  std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to create merge spill: " + out_path.string());
  }

  std::priority_queue<HeapItem, std::vector<HeapItem>, HeapCmp> heap;
  for (std::size_t i = 0; i < cursors.size(); ++i) {
    if (cursors[i].has_value()) {
      heap.push(HeapItem{cursors[i].value(), i});
    }
  }

  std::vector<std::uint64_t> words(n_stripes, 0);
  while (!heap.empty()) {
    const std::uint64_t code = heap.top().code;
    std::fill(words.begin(), words.end(), 0ULL);
    while (!heap.empty() && heap.top().code == code) {
      const HeapItem item = heap.top();
      heap.pop();
      const std::size_t acc = acc_indices[item.src];
      const std::size_t stripe = acc / kAccessionsPerStripe;
      const std::size_t bit = acc % kAccessionsPerStripe;
      words[stripe] |= (1ULL << bit);
      if (auto err = cursors[item.src].advance(); !err.ok()) {
        return err;
      }
      if (cursors[item.src].has_value()) {
        heap.push(HeapItem{cursors[item.src].value(), item.src});
      }
    }
    if (auto err = write_row_file(out, code, words); !err.ok()) {
      return err;
    }
  }
  for (auto& c : cursors) {
    c.close();
  }
  return Error::success();
}

/// Reduce many sorted row spills until count ≤ fd_budget, then return paths.
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

Error process_shard(const fs::path& shard_path, const fs::path& pat_path, const fs::path& map_path,
                    std::size_t n_stripes, std::size_t batch_rows, std::uint64_t& out_patterns,
                    std::uint64_t& out_kmers) {
  const std::size_t record_bytes = sizeof(std::uint64_t) * (1 + n_stripes);
  std::ifstream in(shard_path, std::ios::binary);
  if (!in) {
    return Error::io_error("failed to open shard: " + shard_path.string());
  }
  in.seekg(0, std::ios::end);
  const auto file_bytes = static_cast<std::uint64_t>(in.tellg());
  in.seekg(0, std::ios::beg);
  if (file_bytes % record_bytes != 0) {
    return Error::io_error("corrupt shard size: " + shard_path.string());
  }
  const std::uint64_t nrows = file_bytes / record_bytes;

  std::ofstream pat(pat_path, std::ios::binary | std::ios::trunc);
  std::ofstream map(map_path, std::ios::binary | std::ios::trunc);
  if (!pat || !map) {
    return Error::io_error("failed to open shard outputs");
  }

  std::unordered_map<PatternKey, std::uint32_t, PatternKeyHash> dedup;
  if (nrows > 0) {
    dedup.reserve(static_cast<std::size_t>(std::min<std::uint64_t>(nrows, 1ull << 20)));
  }

  std::vector<char> batch(batch_rows * record_bytes);
  std::vector<std::uint64_t> words(n_stripes);
  std::uint64_t done = 0;
  while (done < nrows) {
    const std::uint64_t take = std::min<std::uint64_t>(batch_rows, nrows - done);
    in.read(batch.data(), static_cast<std::streamsize>(take * record_bytes));
    if (!in) {
      return Error::io_error("failed reading shard rows: " + shard_path.string());
    }
    for (std::uint64_t i = 0; i < take; ++i) {
      const char* rec = batch.data() + i * record_bytes;
      std::uint64_t code = 0;
      std::memcpy(&code, rec, sizeof(code));
      std::memcpy(words.data(), rec + sizeof(code), n_stripes * sizeof(std::uint64_t));

      PatternKey key;
      key.words = words;
      std::uint32_t pid = 0;
      const auto it = dedup.find(key);
      if (it == dedup.end()) {
        if (dedup.size() >= std::numeric_limits<std::uint32_t>::max()) {
          return Error::invalid_argument(
              "pattern count exceeds uint32 in one shard; raise --memory-gb or rebuild with more "
              "partitions");
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
    }
    done += take;
  }

  out_patterns = dedup.size();
  out_kmers = nrows;
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

/// Merge many sorted map shards without holding more than fd_budget FDs.
Error merge_maps_wave(std::vector<fs::path> paths, std::vector<std::uint32_t> offsets,
                      std::size_t fd_budget, const fs::path& work, std::ofstream& out) {
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
      next_offsets.push_back(0);  // already globally adjusted
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
  while (!heap.empty()) {
    const auto item = heap.top();
    heap.pop();
    auto& c = cursors[item.shard];
    out.write(reinterpret_cast<const char*>(&c.ent), sizeof(KmerMapEntry));
    if (!out) {
      return Error::io_error("failed writing k-mer map to matrix");
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

  // Concatenate pattern stores one file at a time (FD-safe for any P).
  const std::size_t batch_rows = opts.batch_rows > 0 ? opts.batch_rows : 100000;
  std::vector<char> copy_buf(std::max<std::size_t>(1u << 20, batch_rows * n_stripes * 8));
  for (std::size_t i = 0; i < num_partitions; ++i) {
    if (pat_counts[i] == 0) {
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
  }

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

  if (auto err = merge_maps_wave(std::move(map_paths), std::move(map_offsets), fd_budget, work, out);
      !err.ok()) {
    return err;
  }

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

  // Header pass: read .kset headers without keeping all FDs (open one at a time).
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

  const std::size_t num_partitions =
      choose_partitions(n_stripes, memory_bytes, threads, batch_rows_req, sum_kmers);
  const std::size_t batch_rows =
      effective_batch_rows(batch_rows_req, num_partitions, n_stripes, memory_bytes);

  {
    std::ostringstream msg;
    msg << "streaming build: accessions=" << n_acc << " stripes=" << n_stripes
        << " sum_kmers=" << sum_kmers << " partitions=" << num_partitions
        << " memory_gb=" << (memory_bytes / (1ull << 30)) << " batch_rows=" << batch_rows
        << " threads=" << threads << " merge_fd_budget=" << fd_budget;
    log_info(msg.str());
  }

  std::vector<ShardWriter> writers(num_partitions);
  for (std::size_t p = 0; p < num_partitions; ++p) {
    const fs::path shard = work / ("shard_" + std::to_string(p) + ".rows");
    if (auto err = writers[p].init(shard, n_stripes, batch_rows); !err.ok()) {
      return err;
    }
  }

  auto partition_row = [&](std::uint64_t code, const std::vector<std::uint64_t>& words) -> Error {
    const std::size_t part = hash_words(words) % num_partitions;
    return writers[part].write_row(code, words);
  };

  std::uint64_t unique_emitted = 0;
  auto count_and_partition = [&](std::uint64_t code,
                                 const std::vector<std::uint64_t>& words) -> Error {
    ++unique_emitted;
    return partition_row(code, words);
  };

  // Direct path: N fits in FD budget — one multiway merge into partitions.
  // Hierarchical path: merge accession groups → row spills → reduce → partitions.
  if (n_acc <= fd_budget) {
    std::vector<PresenceSetCursor> cursors(n_acc);
    for (std::size_t i = 0; i < n_acc; ++i) {
      if (auto err = cursors[i].open(opts.accession_paths[i]); !err.ok()) {
        return err;
      }
    }
    std::priority_queue<HeapItem, std::vector<HeapItem>, HeapCmp> heap;
    for (std::size_t i = 0; i < n_acc; ++i) {
      if (cursors[i].has_value()) {
        heap.push(HeapItem{cursors[i].value(), i});
      }
    }
    std::vector<std::uint64_t> words(n_stripes, 0);
    while (!heap.empty()) {
      const std::uint64_t code = heap.top().code;
      std::fill(words.begin(), words.end(), 0ULL);
      while (!heap.empty() && heap.top().code == code) {
        const HeapItem item = heap.top();
        heap.pop();
        const std::size_t stripe = item.src / kAccessionsPerStripe;
        const std::size_t bit = item.src % kAccessionsPerStripe;
        words[stripe] |= (1ULL << bit);
        if (auto err = cursors[item.src].advance(); !err.ok()) {
          return err;
        }
        if (cursors[item.src].has_value()) {
          heap.push(HeapItem{cursors[item.src].value(), item.src});
        }
      }
      if (auto err = count_and_partition(code, words); !err.ok()) {
        return err;
      }
    }
    for (auto& c : cursors) {
      c.close();
    }
  } else {
    log_info("hierarchical merge: accessions=" + std::to_string(n_acc) +
             " > merge_fd_budget=" + std::to_string(fd_budget));
    const std::size_t group = std::max<std::size_t>(2, fd_budget);
    std::vector<fs::path> spills;
    int spill_seq = 0;
    for (std::size_t begin = 0; begin < n_acc; begin += group) {
      const std::size_t end = std::min(n_acc, begin + group);
      std::vector<std::string> paths;
      std::vector<std::size_t> idxs;
      paths.reserve(end - begin);
      idxs.reserve(end - begin);
      for (std::size_t i = begin; i < end; ++i) {
        paths.push_back(opts.accession_paths[i]);
        idxs.push_back(i);
      }
      const fs::path spill = work / ("acc_merge_" + std::to_string(spill_seq++) + ".rows");
      if (auto err = merge_ksets_to_spill(paths, idxs, n_stripes, spill); !err.ok()) {
        return err;
      }
      spills.push_back(spill);
    }
    if (auto err = reduce_row_spills(spills, n_stripes, fd_budget, work, spill_seq); !err.ok()) {
      return err;
    }
    std::vector<RowCursor> cursors(spills.size());
    for (std::size_t i = 0; i < spills.size(); ++i) {
      if (auto err = cursors[i].open(spills[i], n_stripes); !err.ok()) {
        return err;
      }
    }
    if (auto err = merge_row_cursors(cursors, count_and_partition); !err.ok()) {
      return err;
    }
    for (const auto& p : spills) {
      std::error_code ignore;
      fs::remove(p, ignore);
    }
  }

  for (std::size_t p = 0; p < num_partitions; ++p) {
    if (auto err = writers[p].close(); !err.ok()) {
      return err;
    }
  }

  log_info("merge complete: unique_kmers=" + std::to_string(unique_emitted) +
           "; deduplicating shards");
  std::vector<std::uint64_t> pat_counts(num_partitions, 0);
  std::vector<std::uint64_t> kmer_counts(num_partitions, 0);
  std::vector<Error> errors(num_partitions, Error::success());

  // Concurrent shard workers: each opens ~3 FDs; keep thread count within budget.
  const std::size_t shard_threads = std::min(threads, std::max<std::size_t>(1, fd_budget / 4));

  parallel_for(0, num_partitions, shard_threads, [&](std::size_t p) {
    const fs::path shard = work / ("shard_" + std::to_string(p) + ".rows");
    const fs::path pat = work / ("shard_" + std::to_string(p) + ".pat");
    const fs::path map = work / ("shard_" + std::to_string(p) + ".map");
    std::error_code exists_ec;
    if (!fs::exists(shard, exists_ec) || fs::file_size(shard, exists_ec) == 0) {
      pat_counts[p] = 0;
      kmer_counts[p] = 0;
      std::ofstream(pat, std::ios::binary | std::ios::trunc);
      std::ofstream(map, std::ios::binary | std::ios::trunc);
      return;
    }
    if (auto err = process_shard(shard, pat, map, n_stripes, batch_rows, pat_counts[p],
                                 kmer_counts[p]);
        !err.ok()) {
      errors[p] = err;
    }
    fs::remove(shard, exists_ec);
  });

  for (const Error& err : errors) {
    if (!err.ok()) {
      return err;
    }
  }

  log_info("assembling v2 matrix: patterns=" +
           std::to_string(std::accumulate(pat_counts.begin(), pat_counts.end(), 0ull)) +
           " kmers=" + std::to_string(unique_emitted));

  return assemble_v2(opts, work, n_stripes, num_partitions, pat_counts, kmer_counts, fd_budget);
}

}  // namespace kmat
