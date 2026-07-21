#include "kmat/matrix.hpp"

#include "kmat/fastq.hpp"
#include "kmat/log.hpp"
#include "kmat/presence.hpp"
#include "kmat/runtime.hpp"
#include "kmat/sequence.hpp"
#include "kmat/stripe_build.hpp"
#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace kmat {

namespace {

Error read_header(std::ifstream& in, MatrixHeader& header) {
  in.read(reinterpret_cast<char*>(&header), sizeof(MatrixHeader));
  if (!in) {
    return Error::io_error("failed to read matrix header");
  }
  if (header.magic[0] != 'K' || header.magic[1] != 'M' || header.magic[2] != 'A' ||
      header.magic[3] != 'T') {
    return Error::invalid_argument("invalid matrix magic (expected KMAT)");
  }
  if (header.version != 1 && header.version != 2) {
    return Error::invalid_argument("unsupported matrix version");
  }
  return Error::success();
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

Error inflate_v1_rows(std::vector<MatrixRow> rows, PaMatrix& matrix) {
  matrix.patterns.clear();
  matrix.kmers.clear();
  matrix.kmers.reserve(rows.size());

  std::unordered_map<PatternKey, std::uint32_t, PatternKeyHash> dedup;
  dedup.reserve(rows.size());

  for (MatrixRow& row : rows) {
    PatternKey key;
    key.words = std::move(row.words);
    std::uint32_t pid = 0;
    const auto it = dedup.find(key);
    if (it == dedup.end()) {
      pid = static_cast<std::uint32_t>(matrix.patterns.size());
      dedup.emplace(key, pid);
      matrix.patterns.push_back(key.words);
    } else {
      pid = it->second;
    }
    KmerMapEntry ent;
    ent.kmer_code = row.kmer_code;
    ent.pattern_id = pid;
    matrix.kmers.push_back(ent);
  }

  std::sort(matrix.kmers.begin(), matrix.kmers.end(),
            [](const KmerMapEntry& a, const KmerMapEntry& b) { return a.kmer_code < b.kmer_code; });

  matrix.header.version = 2;
  matrix.header.num_rows = matrix.kmers.size();
  matrix.header.reserved = matrix.patterns.size();
  return Error::success();
}

Error write_v2_from_accession_kmers(const BuildOptions& opts,
                                    std::vector<std::vector<std::uint64_t>> accession_kmers) {
  const std::size_t n_acc = opts.accession_paths.size();
  if (accession_kmers.size() != n_acc) {
    return Error::invalid_argument("accession k-mer sets size mismatch");
  }
  const std::size_t n_stripes = stripe_count_for_accessions(n_acc);

  std::map<std::uint64_t, std::vector<std::uint64_t>> kmer_words;
  for (std::size_t i = 0; i < n_acc; ++i) {
    const std::size_t stripe = i / kAccessionsPerStripe;
    const std::size_t bit = i % kAccessionsPerStripe;
    for (std::uint64_t code : accession_kmers[i]) {
      auto it = kmer_words.find(code);
      if (it == kmer_words.end()) {
        it = kmer_words.emplace(code, std::vector<std::uint64_t>(n_stripes, 0ULL)).first;
      }
      it->second[stripe] |= (1ULL << bit);
    }
  }

  PaMatrix matrix;
  matrix.header.version = 2;
  matrix.header.kmer_size = static_cast<std::uint32_t>(opts.kmer_size);
  matrix.header.num_accessions = static_cast<std::uint32_t>(n_acc);
  matrix.header.num_stripes = static_cast<std::uint32_t>(n_stripes);

  std::unordered_map<PatternKey, std::uint32_t, PatternKeyHash> dedup;
  dedup.reserve(kmer_words.size());
  matrix.kmers.reserve(kmer_words.size());

  for (auto& [code, words] : kmer_words) {
    PatternKey key;
    key.words = words;
    std::uint32_t pid = 0;
    const auto it = dedup.find(key);
    if (it == dedup.end()) {
      pid = static_cast<std::uint32_t>(matrix.patterns.size());
      dedup.emplace(std::move(key), pid);
      matrix.patterns.push_back(std::move(words));
    } else {
      pid = it->second;
    }
    KmerMapEntry ent;
    ent.kmer_code = code;
    ent.pattern_id = pid;
    matrix.kmers.push_back(ent);
  }

  std::sort(matrix.kmers.begin(), matrix.kmers.end(),
            [](const KmerMapEntry& a, const KmerMapEntry& b) { return a.kmer_code < b.kmer_code; });

  matrix.header.num_rows = matrix.kmers.size();
  matrix.header.reserved = matrix.patterns.size();
  return write_matrix(opts.output_path, matrix);
}

Error validate_build_options(const BuildOptions& opts) {
  if (opts.accession_paths.empty()) {
    return Error::invalid_argument("accession list is empty");
  }
  if (opts.kmer_size == 0 || opts.kmer_size > 32) {
    return Error::invalid_argument("k-mer size must be in 1..32");
  }
  if (opts.output_path.empty()) {
    return Error::invalid_argument("output path is required");
  }
  return Error::success();
}

}  // namespace

Error read_matrix_list(const std::string& list_path, std::vector<std::string>& stripe_paths) {
  std::ifstream in(list_path);
  if (!in) {
    return Error::io_error("failed to open matrix list: " + list_path);
  }
  stripe_paths.clear();
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty() || line.front() == '#') {
      continue;
    }
    stripe_paths.push_back(line);
  }
  return Error::success();
}

Error read_matrix(const std::string& path, PaMatrix& matrix) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return Error::io_error("failed to open matrix: " + path);
  }

  matrix = PaMatrix{};
  if (auto err = read_header(in, matrix.header); !err.ok()) {
    return err;
  }

  if (matrix.header.version == 1) {
    std::vector<MatrixRow> rows(static_cast<std::size_t>(matrix.header.num_rows));
    for (std::size_t i = 0; i < rows.size(); ++i) {
      MatrixRow& row = rows[i];
      row.words.resize(matrix.header.num_stripes);
      in.read(reinterpret_cast<char*>(&row.kmer_code), sizeof(std::uint64_t));
      if (!in) {
        return Error::io_error("unexpected EOF reading kmer code");
      }
      in.read(reinterpret_cast<char*>(row.words.data()),
              static_cast<std::streamsize>(matrix.header.num_stripes * sizeof(std::uint64_t)));
      if (!in) {
        return Error::io_error("unexpected EOF reading presence words");
      }
    }
    return inflate_v1_rows(std::move(rows), matrix);
  }

  // version 2
  const std::uint64_t num_patterns = matrix.header.reserved;
  const std::uint64_t num_kmers = matrix.header.num_rows;
  matrix.patterns.resize(static_cast<std::size_t>(num_patterns));
  for (std::size_t p = 0; p < matrix.patterns.size(); ++p) {
    matrix.patterns[p].resize(matrix.header.num_stripes);
    in.read(reinterpret_cast<char*>(matrix.patterns[p].data()),
            static_cast<std::streamsize>(matrix.header.num_stripes * sizeof(std::uint64_t)));
    if (!in) {
      return Error::io_error("unexpected EOF reading pattern store");
    }
  }

  matrix.kmers.resize(static_cast<std::size_t>(num_kmers));
  in.read(reinterpret_cast<char*>(matrix.kmers.data()),
          static_cast<std::streamsize>(num_kmers * sizeof(KmerMapEntry)));
  if (!in) {
    return Error::io_error("unexpected EOF reading k-mer map");
  }

  return Error::success();
}

Error write_matrix(const std::string& path, const PaMatrix& matrix) {
  if (matrix.header.version != 2) {
    return Error::invalid_argument("write_matrix currently emits version 2 only");
  }
  if (matrix.patterns.size() != matrix.header.reserved) {
    return Error::invalid_argument("pattern count does not match header.reserved");
  }
  if (matrix.kmers.size() != matrix.header.num_rows) {
    return Error::invalid_argument("k-mer count does not match header.num_rows");
  }

  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to open matrix for writing: " + path);
  }

  const std::size_t buf_sz = std::max<std::size_t>(1u << 16, runtime_config().io_buffer_bytes);
  std::vector<char> buf(buf_sz);
  out.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));

  out.write(reinterpret_cast<const char*>(&matrix.header), sizeof(MatrixHeader));
  for (const auto& words : matrix.patterns) {
    if (words.size() != matrix.header.num_stripes) {
      return Error::invalid_argument("pattern stripe width mismatch");
    }
    out.write(reinterpret_cast<const char*>(words.data()),
              static_cast<std::streamsize>(words.size() * sizeof(std::uint64_t)));
  }
  out.write(reinterpret_cast<const char*>(matrix.kmers.data()),
            static_cast<std::streamsize>(matrix.kmers.size() * sizeof(KmerMapEntry)));

  if (!out) {
    return Error::io_error("failed while writing matrix: " + path);
  }
  return Error::success();
}

Error load_matrix_from_list(const std::string& matrix_list_path, std::size_t num_accessions,
                            PaMatrix& matrix) {
  std::vector<std::string> stripes;
  if (auto err = read_matrix_list(matrix_list_path, stripes); !err.ok()) {
    return err;
  }
  if (stripes.empty()) {
    return Error::invalid_argument("matrix list is empty");
  }

  const std::size_t expected_stripes = stripe_count_for_accessions(num_accessions);
  if (stripes.size() != expected_stripes) {
    return Error::invalid_argument("stripe count mismatch in matrix list");
  }

  // Legacy multi-file stripes are v1 dense single-word files. Merge then inflate.
  std::ifstream first(stripes.front(), std::ios::binary);
  if (!first) {
    return Error::io_error("failed to open stripe: " + stripes.front());
  }
  MatrixHeader header{};
  if (auto err = read_header(first, header); !err.ok()) {
    return err;
  }
  if (header.version != 1) {
    return Error::invalid_argument("matrix list loading supports version 1 stripe files only");
  }
  if (header.num_accessions != num_accessions) {
    return Error::invalid_argument("matrix accession count does not match accession list");
  }

  std::vector<MatrixRow> rows(static_cast<std::size_t>(header.num_rows));
  for (std::size_t i = 0; i < rows.size(); ++i) {
    rows[i].words.assign(expected_stripes, 0ULL);
    first.read(reinterpret_cast<char*>(&rows[i].kmer_code), sizeof(std::uint64_t));
    std::uint64_t w0 = 0;
    first.read(reinterpret_cast<char*>(&w0), sizeof(std::uint64_t));
    if (!first) {
      return Error::io_error("EOF reading first stripe");
    }
    rows[i].words[0] = w0;
  }
  first.close();

  for (std::size_t s = 1; s < stripes.size(); ++s) {
    std::ifstream in(stripes[s], std::ios::binary);
    if (!in) {
      return Error::io_error("failed to open stripe: " + stripes[s]);
    }
    MatrixHeader h{};
    if (auto err = read_header(in, h); !err.ok()) {
      return err;
    }
    if (h.num_rows != header.num_rows || h.kmer_size != header.kmer_size) {
      return Error::invalid_argument("stripe row metadata mismatch");
    }
    for (std::size_t r = 0; r < rows.size(); ++r) {
      std::uint64_t code = 0;
      std::uint64_t w = 0;
      in.read(reinterpret_cast<char*>(&code), sizeof(std::uint64_t));
      in.read(reinterpret_cast<char*>(&w), sizeof(std::uint64_t));
      if (!in) {
        return Error::io_error("EOF reading stripe");
      }
      if (code != rows[r].kmer_code) {
        return Error::invalid_argument("k-mer row order mismatch across stripes");
      }
      rows[r].words[s] = w;
    }
  }

  matrix.header = header;
  matrix.header.num_stripes = static_cast<std::uint32_t>(expected_stripes);
  return inflate_v1_rows(std::move(rows), matrix);
}

Error build_matrix_from_sequences(const BuildOptions& opts) {
  if (auto err = validate_build_options(opts); !err.ok()) {
    return err;
  }

  log_warn(
      "building from sequence files loads k-mers in RAM; production panels should use "
      "`kmat count` → `.kset` then memory-bounded `kmat build`");

  const std::size_t n_acc = opts.accession_paths.size();
  std::vector<std::vector<std::uint64_t>> accession_kmers(n_acc);
  std::vector<Error> errors(n_acc, Error::success());
  const std::size_t threads = opts.num_threads > 0 ? opts.num_threads
                                                   : effective_threads(runtime_config());

  parallel_for(0, n_acc, threads, [&](std::size_t i) {
    std::unordered_set<std::uint64_t> uniq;
    std::vector<FastaRecord> records;
    if (auto err = read_sequence_file(opts.accession_paths[i], records); !err.ok()) {
      errors[i] = err;
      return;
    }
    for (const FastaRecord& rec : records) {
      if (auto err = kmer_set_from_sequence(rec.sequence, opts.kmer_size, uniq); !err.ok()) {
        errors[i] = err;
        return;
      }
    }
    accession_kmers[i].assign(uniq.begin(), uniq.end());
  });

  for (const Error& err : errors) {
    if (!err.ok()) {
      return err;
    }
  }

  return write_v2_from_accession_kmers(opts, std::move(accession_kmers));
}

Error build_matrix_from_presence_sets(const BuildOptions& opts) {
  if (auto err = validate_build_options(opts); !err.ok()) {
    return err;
  }
  // Production path: master → stripe create/fill → v2 compress (few large files).
  return build_matrix_from_presence_sets_staged(opts);
}

Error build_matrix_from_accessions(const BuildOptions& opts) {
  if (opts.accession_paths.empty()) {
    return Error::invalid_argument("accession list is empty");
  }
  std::size_t n_kset = 0;
  for (const std::string& path : opts.accession_paths) {
    if (path_looks_presence_set(path)) {
      ++n_kset;
    }
  }
  if (n_kset == opts.accession_paths.size()) {
    return build_matrix_from_presence_sets(opts);
  }
  if (n_kset == 0) {
    return build_matrix_from_sequences(opts);
  }
  return Error::invalid_argument(
      "accession list mixes .kset presence sets with sequence files; use one input type");
}

Error build_matrix_from_fastas(const BuildOptions& opts) {
  return build_matrix_from_sequences(opts);
}

bool get_presence_bit(const std::vector<std::uint64_t>& words, std::size_t accession_index) {
  const std::size_t stripe = accession_index / kAccessionsPerStripe;
  const std::size_t bit = accession_index % kAccessionsPerStripe;
  if (stripe >= words.size()) {
    return false;
  }
  return (words[stripe] & (1ULL << bit)) != 0ULL;
}

bool get_presence_bit(const MatrixRow& row, std::size_t accession_index) {
  return get_presence_bit(row.words, accession_index);
}

std::string presence_bitstring(const std::vector<std::uint64_t>& words, std::size_t num_accessions) {
  std::string bits;
  bits.reserve(num_accessions);
  for (std::size_t i = 0; i < num_accessions; ++i) {
    bits.push_back(get_presence_bit(words, i) ? '1' : '0');
  }
  return bits;
}

std::string presence_bitstring(const MatrixRow& row, std::size_t num_accessions) {
  return presence_bitstring(row.words, num_accessions);
}

std::size_t count_carriers(const std::vector<std::uint64_t>& words, std::size_t num_accessions) {
  std::size_t n = 0;
  for (std::size_t i = 0; i < num_accessions; ++i) {
    n += get_presence_bit(words, i) ? 1 : 0;
  }
  return n;
}

std::size_t count_carriers(const MatrixRow& row, std::size_t num_accessions) {
  return count_carriers(row.words, num_accessions);
}

std::optional<std::uint32_t> find_pattern_id(const PaMatrix& matrix, std::uint64_t kmer_code) {
  const auto it =
      std::lower_bound(matrix.kmers.begin(), matrix.kmers.end(), kmer_code,
                       [](const KmerMapEntry& e, std::uint64_t code) { return e.kmer_code < code; });
  if (it == matrix.kmers.end() || it->kmer_code != kmer_code) {
    return std::nullopt;
  }
  return it->pattern_id;
}

std::vector<std::uint64_t> kmers_for_pattern(const PaMatrix& matrix, std::uint32_t pattern_id) {
  std::vector<std::uint64_t> out;
  for (const KmerMapEntry& e : matrix.kmers) {
    if (e.pattern_id == pattern_id) {
      out.push_back(e.kmer_code);
    }
  }
  return out;
}

std::vector<std::vector<std::uint64_t>> pattern_kmer_index(const PaMatrix& matrix) {
  std::vector<std::vector<std::uint64_t>> index(matrix.patterns.size());
  for (const KmerMapEntry& e : matrix.kmers) {
    if (e.pattern_id < index.size()) {
      index[e.pattern_id].push_back(e.kmer_code);
    }
  }
  return index;
}

}  // namespace kmat
