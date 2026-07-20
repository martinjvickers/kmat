#include "kmat/build_stream.hpp"

#include "kmat/log.hpp"
#include "kmat/matrix_layout.hpp"
#include "kmat/presence.hpp"
#include "kmat/runtime.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
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

std::size_t choose_partitions(std::uint64_t sum_kmers, std::size_t n_stripes,
                              std::size_t memory_bytes, std::size_t num_threads) {
  // Conservative: treat sum_kmers as an upper bound on unique U; size each shard so
  // pattern table + row buffer fits in (memory / threads).
  const std::size_t threads = std::max<std::size_t>(1, num_threads);
  const std::size_t budget = std::max<std::size_t>(1ull << 20, memory_bytes / threads);
  // ~24B overhead + code + words + pattern key copy during dedup
  const std::size_t bytes_per_unique =
      64 + (1 + n_stripes) * sizeof(std::uint64_t) + n_stripes * sizeof(std::uint64_t);
  const std::uint64_t per_shard_cap =
      std::max<std::uint64_t>(1, static_cast<std::uint64_t>(budget / bytes_per_unique));
  std::uint64_t p = (sum_kmers + per_shard_cap - 1) / per_shard_cap;
  p = std::max<std::uint64_t>(p, threads);
  p = std::min<std::uint64_t>(p, 4096);
  p = std::max<std::uint64_t>(p, 1);
  return static_cast<std::size_t>(p);
}

struct HeapItem {
  std::uint64_t code{0};
  std::size_t acc{0};
};

struct HeapCmp {
  bool operator()(const HeapItem& a, const HeapItem& b) const {
    if (a.code != b.code) {
      return a.code > b.code;  // min-heap by code
    }
    return a.acc > b.acc;
  }
};

struct ShardWriter {
  std::ofstream out;
  std::vector<char> buf;
  std::size_t buffered_rows{0};
  std::size_t batch_rows{0};
  std::size_t record_bytes{0};
  std::uint64_t total_rows{0};

  Error open(const fs::path& path, std::size_t n_stripes, std::size_t batch) {
    record_bytes = sizeof(std::uint64_t) * (1 + n_stripes);
    batch_rows = std::max<std::size_t>(1, batch);
    buf.reserve(batch_rows * record_bytes);
    out.open(path, std::ios::binary | std::ios::trunc);
    if (!out) {
      return Error::io_error("failed to open shard for writing: " + path.string());
    }
    return Error::success();
  }

  Error flush() {
    if (buf.empty()) {
      return Error::success();
    }
    out.write(buf.data(), static_cast<std::streamsize>(buf.size()));
    if (!out) {
      return Error::io_error("failed writing shard rows");
    }
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

  Error close() {
    if (auto err = flush(); !err.ok()) {
      return err;
    }
    out.close();
    return Error::success();
  }
};

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

Error assemble_v2(const BuildOptions& opts, const fs::path& work, std::size_t n_stripes,
                  std::size_t num_partitions, const std::vector<std::uint64_t>& pat_counts,
                  const std::vector<std::uint64_t>& kmer_counts) {
  std::uint64_t total_patterns = 0;
  std::uint64_t total_kmers = 0;
  std::vector<std::uint32_t> offsets(num_partitions, 0);
  for (std::size_t i = 0; i < num_partitions; ++i) {
    offsets[i] = static_cast<std::uint32_t>(total_patterns);
    total_patterns += pat_counts[i];
    total_kmers += kmer_counts[i];
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

  // Concatenate pattern stores in shard order.
  std::vector<char> copy_buf(std::max<std::size_t>(1u << 20, opts.batch_rows * n_stripes * 8));
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

  // Multiway merge map shards (each sorted by code) with global pattern ids.
  std::vector<MapCursor> cursors(num_partitions);
  std::priority_queue<MapHeapItem, std::vector<MapHeapItem>, MapHeapCmp> heap;
  for (std::size_t i = 0; i < num_partitions; ++i) {
    if (kmer_counts[i] == 0) {
      continue;
    }
    const fs::path map = work / ("shard_" + std::to_string(i) + ".map");
    if (auto err = cursors[i].open(map, offsets[i]); !err.ok()) {
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
  const std::size_t batch_rows = opts.batch_rows > 0 ? opts.batch_rows : 100000;

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

  std::vector<PresenceSetCursor> cursors(n_acc);
  std::uint64_t sum_kmers = 0;
  for (std::size_t i = 0; i < n_acc; ++i) {
    if (auto err = cursors[i].open(opts.accession_paths[i]); !err.ok()) {
      return err;
    }
    if (cursors[i].header().kmer_size != opts.kmer_size) {
      return Error::invalid_argument("presence set k-mer size mismatch: " +
                                     opts.accession_paths[i]);
    }
    sum_kmers += cursors[i].header().num_kmers;
  }

  const std::size_t num_partitions =
      choose_partitions(sum_kmers, n_stripes, memory_bytes, threads);

  {
    std::ostringstream msg;
    msg << "streaming build: accessions=" << n_acc << " stripes=" << n_stripes
        << " sum_kmers=" << sum_kmers << " partitions=" << num_partitions
        << " memory_gb=" << (memory_bytes / (1ull << 30)) << " batch_rows=" << batch_rows
        << " threads=" << threads;
    log_info(msg.str());
  }

  std::vector<ShardWriter> writers(num_partitions);
  for (std::size_t p = 0; p < num_partitions; ++p) {
    const fs::path shard = work / ("shard_" + std::to_string(p) + ".rows");
    if (auto err = writers[p].open(shard, n_stripes, batch_rows); !err.ok()) {
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
  std::uint64_t unique_emitted = 0;
  while (!heap.empty()) {
    const std::uint64_t code = heap.top().code;
    std::fill(words.begin(), words.end(), 0ULL);
    while (!heap.empty() && heap.top().code == code) {
      const HeapItem item = heap.top();
      heap.pop();
      const std::size_t stripe = item.acc / kAccessionsPerStripe;
      const std::size_t bit = item.acc % kAccessionsPerStripe;
      words[stripe] |= (1ULL << bit);
      if (auto err = cursors[item.acc].advance(); !err.ok()) {
        return err;
      }
      if (cursors[item.acc].has_value()) {
        heap.push(HeapItem{cursors[item.acc].value(), item.acc});
      }
    }

    const std::size_t part = hash_words(words) % num_partitions;
    if (auto err = writers[part].write_row(code, words); !err.ok()) {
      return err;
    }
    ++unique_emitted;
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

  parallel_for(0, num_partitions, threads, [&](std::size_t p) {
    const fs::path shard = work / ("shard_" + std::to_string(p) + ".rows");
    const fs::path pat = work / ("shard_" + std::to_string(p) + ".pat");
    const fs::path map = work / ("shard_" + std::to_string(p) + ".map");
    std::error_code exists_ec;
    if (!fs::exists(shard, exists_ec) || fs::file_size(shard, exists_ec) == 0) {
      pat_counts[p] = 0;
      kmer_counts[p] = 0;
      // Touch empty outputs for uniform assemble.
      std::ofstream(pat, std::ios::binary | std::ios::trunc);
      std::ofstream(map, std::ios::binary | std::ios::trunc);
      return;
    }
    if (auto err = process_shard(shard, pat, map, n_stripes, batch_rows, pat_counts[p],
                                 kmer_counts[p]);
        !err.ok()) {
      errors[p] = err;
    }
    // Free shard row file early.
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

  return assemble_v2(opts, work, n_stripes, num_partitions, pat_counts, kmer_counts);
}

}  // namespace kmat
