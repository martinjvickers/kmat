#include "kmat/stripe_build.hpp"

#include "kmat/log.hpp"
#include "kmat/matrix_layout.hpp"
#include "kmat/presence.hpp"
#include "kmat/runtime.hpp"
#include "kmat/universe.hpp"

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <unistd.h>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace kmat {

namespace {

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

std::size_t default_memory_bytes() {
  const auto cfg = runtime_config();
  if (cfg.profile == RuntimeProfile::Hpc) {
    return 64ull << 30;
  }
  return 8ull << 30;
}

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

Error write_v1_stripe_header(std::ofstream& out, std::size_t kmer_size, std::size_t num_accessions,
                             std::uint64_t num_rows) {
  MatrixHeader header{};
  header.magic[0] = 'K';
  header.magic[1] = 'M';
  header.magic[2] = 'A';
  header.magic[3] = 'T';
  header.version = 1;
  header.kmer_size = static_cast<std::uint32_t>(kmer_size);
  header.num_accessions = static_cast<std::uint32_t>(num_accessions);
  header.num_stripes = 1;  // one word per row in this file
  header.num_rows = num_rows;
  header.reserved = 0;
  out.write(reinterpret_cast<const char*>(&header), sizeof(header));
  if (!out) {
    return Error::io_error("failed writing stripe header");
  }
  return Error::success();
}

}  // namespace

Error create_blank_stripe(const CreateStripeOptions& opts) {
  if (opts.universe_path.empty() || opts.output_path.empty()) {
    return Error::invalid_argument("universe and output paths are required");
  }
  if (opts.num_accessions == 0) {
    return Error::invalid_argument("num_accessions must be > 0");
  }
  const std::size_t n_stripes = stripe_count_for_accessions(opts.num_accessions);
  if (opts.stripe_index >= n_stripes) {
    return Error::invalid_argument("stripe_index out of range");
  }

  UniverseCursor uni;
  if (auto err = uni.open(opts.universe_path); !err.ok()) {
    return err;
  }
  if (uni.header().kmer_size != opts.kmer_size) {
    return Error::invalid_argument("universe k-mer size mismatch");
  }

  const std::size_t batch = opts.batch_rows > 0 ? opts.batch_rows : kDefaultBatchRows;
  constexpr std::size_t kRecordBytes = 16;  // code + word

  std::ofstream out(opts.output_path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to open stripe for writing: " + opts.output_path);
  }
  if (auto err = write_v1_stripe_header(out, opts.kmer_size, opts.num_accessions,
                                        uni.header().num_kmers);
      !err.ok()) {
    return err;
  }

  log_info("create-stripe: stripe=" + std::to_string(opts.stripe_index) + " rows=" +
           std::to_string(uni.header().num_kmers) + " batch=" + std::to_string(batch));

  std::vector<char> buf(batch * kRecordBytes);
  std::uint64_t written = 0;
  const auto t0 = std::chrono::steady_clock::now();
  std::size_t buffered = 0;

  auto flush = [&]() -> Error {
    if (buffered == 0) {
      return Error::success();
    }
    out.write(buf.data(), static_cast<std::streamsize>(buffered * kRecordBytes));
    if (!out) {
      return Error::io_error("failed writing stripe rows");
    }
    buffered = 0;
    return Error::success();
  };

  while (uni.has_value()) {
    const std::size_t off = buffered * kRecordBytes;
    const std::uint64_t code = uni.value();
    const std::uint64_t word = 0;
    std::memcpy(buf.data() + off, &code, sizeof(code));
    std::memcpy(buf.data() + off + sizeof(code), &word, sizeof(word));
    ++buffered;
    ++written;
    if (buffered >= batch) {
      if (auto err = flush(); !err.ok()) {
        return err;
      }
      const auto elapsed =
          std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0)
              .count();
      if (elapsed > 0 && (written % (batch * 50) == 0)) {
        log_info("create-stripe: rows=" + std::to_string(written) + "/" +
                 std::to_string(uni.header().num_kmers) + " elapsed=" + std::to_string(elapsed) +
                 "s");
      }
    }
    if (auto err = uni.advance(); !err.ok()) {
      return err;
    }
  }
  if (auto err = flush(); !err.ok()) {
    return err;
  }
  uni.close();
  log_info("create-stripe: done rows=" + std::to_string(written) + " path=" + opts.output_path);
  return Error::success();
}

Error fill_stripe_column(const FillStripeOptions& opts) {
  if (opts.stripe_path.empty() || opts.kset_path.empty()) {
    return Error::invalid_argument("stripe and kset paths are required");
  }
  if (opts.local_index >= kAccessionsPerStripe) {
    return Error::invalid_argument("local_index must be in 0..63");
  }

  PresenceSetCursor kset;
  if (auto err = kset.open(opts.kset_path); !err.ok()) {
    return err;
  }

  std::fstream io(opts.stripe_path, std::ios::in | std::ios::out | std::ios::binary);
  if (!io) {
    return Error::io_error("failed to open stripe for fill: " + opts.stripe_path);
  }
  MatrixHeader header{};
  io.read(reinterpret_cast<char*>(&header), sizeof(header));
  if (!io) {
    return Error::io_error("failed reading stripe header");
  }
  if (header.magic[0] != 'K' || header.magic[1] != 'M' || header.magic[2] != 'A' ||
      header.magic[3] != 'T' || header.version != 1) {
    return Error::invalid_argument("fill expects a v1 dense stripe file");
  }
  if (header.num_stripes != 1) {
    return Error::invalid_argument("fill expects single-word stripe files");
  }
  if (kset.header().kmer_size != header.kmer_size) {
    return Error::invalid_argument("kset k-mer size mismatch vs stripe");
  }

  const std::size_t batch = opts.batch_rows > 0 ? opts.batch_rows : kDefaultBatchRows;
  constexpr std::size_t kRecordBytes = 16;
  std::vector<char> buf(batch * kRecordBytes);
  const std::uint64_t bit = (1ULL << opts.local_index);

  log_info("fill: local_index=" + std::to_string(opts.local_index) + " rows=" +
           std::to_string(header.num_rows) + " batch=" + std::to_string(batch));

  std::uint64_t index = 0;
  const auto t0 = std::chrono::steady_clock::now();
  while (index < header.num_rows) {
    const std::uint64_t take =
        std::min<std::uint64_t>(batch, header.num_rows - index);
    const std::uint64_t file_offset = sizeof(MatrixHeader) + index * kRecordBytes;
    io.seekg(static_cast<std::streamoff>(file_offset));
    io.read(buf.data(), static_cast<std::streamsize>(take * kRecordBytes));
    if (!io) {
      return Error::io_error("failed reading stripe batch");
    }

    for (std::uint64_t i = 0; i < take; ++i) {
      std::uint64_t code = 0;
      std::memcpy(&code, buf.data() + i * kRecordBytes, sizeof(code));
      // Advance kset to >= code
      while (kset.has_value() && kset.value() < code) {
        if (auto err = kset.advance(); !err.ok()) {
          return err;
        }
      }
      if (kset.has_value() && kset.value() == code) {
        std::uint64_t word = 0;
        std::memcpy(&word, buf.data() + i * kRecordBytes + sizeof(code), sizeof(word));
        word |= bit;
        std::memcpy(buf.data() + i * kRecordBytes + sizeof(code), &word, sizeof(word));
      }
    }

    io.seekp(static_cast<std::streamoff>(file_offset));
    io.write(buf.data(), static_cast<std::streamsize>(take * kRecordBytes));
    if (!io) {
      return Error::io_error("failed writing stripe batch");
    }
    index += take;

    const auto elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0)
            .count();
    if (elapsed >= 30 && (index % (batch * 20) == 0 || index == header.num_rows)) {
      log_info("fill: rows=" + std::to_string(index) + "/" + std::to_string(header.num_rows) +
               " elapsed=" + std::to_string(elapsed) + "s");
    }
  }

  kset.close();
  log_info("fill: done local_index=" + std::to_string(opts.local_index));
  return Error::success();
}

Error compress_stripes_to_v2(const CompressOptions& opts) {
  if (opts.stripe_paths.empty() || opts.output_path.empty()) {
    return Error::invalid_argument("stripe list and output are required");
  }
  if (opts.num_accessions == 0) {
    return Error::invalid_argument("num_accessions must be > 0");
  }
  const std::size_t n_stripes = stripe_count_for_accessions(opts.num_accessions);
  if (opts.stripe_paths.size() != n_stripes) {
    return Error::invalid_argument("stripe count mismatch vs num_accessions");
  }

  std::vector<std::ifstream> ins(n_stripes);
  MatrixHeader first_hdr{};
  for (std::size_t s = 0; s < n_stripes; ++s) {
    ins[s].open(opts.stripe_paths[s], std::ios::binary);
    if (!ins[s]) {
      return Error::io_error("failed to open stripe: " + opts.stripe_paths[s]);
    }
    MatrixHeader h{};
    ins[s].read(reinterpret_cast<char*>(&h), sizeof(h));
    if (!ins[s] || h.version != 1 || h.num_stripes != 1) {
      return Error::invalid_argument("compress expects v1 single-word stripes");
    }
    if (s == 0) {
      first_hdr = h;
    } else if (h.num_rows != first_hdr.num_rows || h.kmer_size != first_hdr.kmer_size) {
      return Error::invalid_argument("stripe metadata mismatch");
    }
    if (h.kmer_size != opts.kmer_size) {
      return Error::invalid_argument("stripe k-mer size mismatch");
    }
  }

  const std::size_t batch = opts.batch_rows > 0 ? opts.batch_rows : kDefaultBatchRows;
  const std::size_t memory_bytes =
      opts.memory_bytes > 0 ? opts.memory_bytes : default_memory_bytes();

  log_info("compress: rows=" + std::to_string(first_hdr.num_rows) + " stripes=" +
           std::to_string(n_stripes) + " memory_gb=" + std::to_string(memory_bytes >> 30));

  const fs::path work =
      fs::path(scratch_root(opts.tmpdir)) / ("kmat_compress_" + std::to_string(::getpid()));
  std::error_code ec;
  fs::create_directories(work, ec);
  struct Cleaner {
    fs::path path;
    ~Cleaner() {
      std::error_code ignore;
      fs::remove_all(path, ignore);
    }
  } cleaner{work};

  const fs::path pat_path = work / "patterns.bin";
  const fs::path map_path = work / "map.bin";
  std::ofstream pat_out(pat_path, std::ios::binary | std::ios::trunc);
  std::ofstream map_out(map_path, std::ios::binary | std::ios::trunc);
  if (!pat_out || !map_out) {
    return Error::io_error("failed to open compress spills");
  }

  std::unordered_map<PatternKey, std::uint32_t, PatternKeyHash> dedup;
  dedup.reserve(1u << 20);

  std::vector<std::uint64_t> words(n_stripes);
  std::vector<std::uint64_t> codes(batch);
  std::vector<std::uint64_t> batch_words(static_cast<std::size_t>(batch) * n_stripes);
  std::uint64_t done = 0;
  const auto t0 = std::chrono::steady_clock::now();

  while (done < first_hdr.num_rows) {
    const std::uint64_t take =
        std::min<std::uint64_t>(batch, first_hdr.num_rows - done);
    for (std::uint64_t i = 0; i < take; ++i) {
      std::uint64_t code0 = 0;
      for (std::size_t s = 0; s < n_stripes; ++s) {
        std::uint64_t code = 0;
        std::uint64_t w = 0;
        ins[s].read(reinterpret_cast<char*>(&code), sizeof(code));
        ins[s].read(reinterpret_cast<char*>(&w), sizeof(w));
        if (!ins[s]) {
          return Error::io_error("EOF reading stripe during compress");
        }
        if (s == 0) {
          code0 = code;
          codes[static_cast<std::size_t>(i)] = code;
        } else if (code != code0) {
          return Error::invalid_argument("k-mer order mismatch across stripes");
        }
        batch_words[static_cast<std::size_t>(i) * n_stripes + s] = w;
      }
    }

    for (std::uint64_t i = 0; i < take; ++i) {
      for (std::size_t s = 0; s < n_stripes; ++s) {
        words[s] = batch_words[static_cast<std::size_t>(i) * n_stripes + s];
      }
      PatternKey key;
      key.words = words;
      std::uint32_t pid = 0;
      const auto it = dedup.find(key);
      if (it == dedup.end()) {
        if (dedup.size() >= std::numeric_limits<std::uint32_t>::max()) {
          return Error::invalid_argument("pattern count exceeds uint32");
        }
        pid = static_cast<std::uint32_t>(dedup.size());
        dedup.emplace(std::move(key), pid);
        pat_out.write(reinterpret_cast<const char*>(words.data()),
                      static_cast<std::streamsize>(n_stripes * sizeof(std::uint64_t)));
        if (!pat_out) {
          return Error::io_error("failed writing pattern store");
        }
      } else {
        pid = it->second;
      }
      KmerMapEntry ent;
      ent.kmer_code = codes[static_cast<std::size_t>(i)];
      ent.pattern_id = pid;
      ent.pad = 0;
      map_out.write(reinterpret_cast<const char*>(&ent), sizeof(ent));
      if (!map_out) {
        return Error::io_error("failed writing k-mer map");
      }
    }

    done += take;
    const auto elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0)
            .count();
    if (elapsed >= 30 || done == first_hdr.num_rows) {
      log_info("compress: rows=" + std::to_string(done) + "/" +
               std::to_string(first_hdr.num_rows) + " patterns=" + std::to_string(dedup.size()) +
               " elapsed=" + std::to_string(elapsed) + "s");
    }
  }

  pat_out.close();
  map_out.close();

  MatrixHeader out_hdr{};
  out_hdr.magic[0] = 'K';
  out_hdr.magic[1] = 'M';
  out_hdr.magic[2] = 'A';
  out_hdr.magic[3] = 'T';
  out_hdr.version = 2;
  out_hdr.kmer_size = static_cast<std::uint32_t>(opts.kmer_size);
  out_hdr.num_accessions = static_cast<std::uint32_t>(opts.num_accessions);
  out_hdr.num_stripes = static_cast<std::uint32_t>(n_stripes);
  out_hdr.num_rows = first_hdr.num_rows;
  out_hdr.reserved = dedup.size();

  std::ofstream out(opts.output_path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to open v2 output: " + opts.output_path);
  }
  out.write(reinterpret_cast<const char*>(&out_hdr), sizeof(out_hdr));

  // Concat pattern store then map
  auto copy_file = [&](const fs::path& p) -> Error {
    std::ifstream in(p, std::ios::binary);
    if (!in) {
      return Error::io_error("failed to reopen spill: " + p.string());
    }
    std::vector<char> copy_buf(1u << 20);
    while (in) {
      in.read(copy_buf.data(), static_cast<std::streamsize>(copy_buf.size()));
      const auto n = in.gcount();
      if (n > 0) {
        out.write(copy_buf.data(), n);
        if (!out) {
          return Error::io_error("failed writing v2 body");
        }
      }
    }
    return Error::success();
  };
  if (auto err = copy_file(pat_path); !err.ok()) {
    return err;
  }
  if (auto err = copy_file(map_path); !err.ok()) {
    return err;
  }

  log_info("compress: done patterns=" + std::to_string(dedup.size()) + " kmers=" +
           std::to_string(first_hdr.num_rows) + " path=" + opts.output_path);
  return Error::success();
}

Error build_matrix_from_presence_sets_staged(const BuildOptions& opts) {
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
  const std::size_t batch =
      opts.batch_rows > 0 ? opts.batch_rows : kDefaultBatchRows;
  const std::size_t memory_bytes =
      opts.memory_bytes > 0 ? opts.memory_bytes : default_memory_bytes();

  const fs::path work =
      fs::path(scratch_root(opts.tmpdir)) / ("kmat_staged_" + std::to_string(::getpid()));
  std::error_code ec;
  fs::create_directories(work, ec);
  if (ec) {
    return Error::io_error("failed to create staged build work dir: " + work.string());
  }
  struct Cleaner {
    fs::path path;
    ~Cleaner() {
      std::error_code ignore;
      fs::remove_all(path, ignore);
    }
  } cleaner{work};

  log_info("staged build: accessions=" + std::to_string(n_acc) + " stripes=" +
           std::to_string(n_stripes) + " batch_rows=" + std::to_string(batch) +
           " memory_gb=" + std::to_string(memory_bytes >> 30));

  const fs::path kuniv = work / "panel.kuniv";
  if (auto err = build_universe_from_presence_sets(opts.accession_paths, opts.kmer_size,
                                                   kuniv.string(), work.string(), 32);
      !err.ok()) {
    return err;
  }

  std::vector<std::string> stripe_paths;
  stripe_paths.reserve(n_stripes);
  for (std::size_t s = 0; s < n_stripes; ++s) {
    const fs::path stripe = work / ("panel." + std::to_string(s) + ".bin");
    CreateStripeOptions copts;
    copts.universe_path = kuniv.string();
    copts.output_path = stripe.string();
    copts.kmer_size = opts.kmer_size;
    copts.num_accessions = n_acc;
    copts.stripe_index = s;
    copts.batch_rows = batch;
    if (auto err = create_blank_stripe(copts); !err.ok()) {
      return err;
    }
    stripe_paths.push_back(stripe.string());
  }

  for (std::size_t g = 0; g < n_acc; ++g) {
    const std::size_t s = g / kAccessionsPerStripe;
    const std::size_t local = g % kAccessionsPerStripe;
    FillStripeOptions fopts;
    fopts.stripe_path = stripe_paths[s];
    fopts.kset_path = opts.accession_paths[g];
    fopts.local_index = local;
    fopts.batch_rows = batch;
    if (auto err = fill_stripe_column(fopts); !err.ok()) {
      return err;
    }
  }

  CompressOptions zopts;
  zopts.stripe_paths = stripe_paths;
  zopts.num_accessions = n_acc;
  zopts.kmer_size = opts.kmer_size;
  zopts.output_path = opts.output_path;
  zopts.memory_bytes = memory_bytes;
  zopts.batch_rows = batch;
  zopts.tmpdir = work.string();
  return compress_stripes_to_v2(zopts);
}

}  // namespace kmat
