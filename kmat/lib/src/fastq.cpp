#include "kmat/fastq.hpp"

#include "kmat/fasta.hpp"

#include <cctype>
#include <fstream>
#include <functional>
#include <string>
#include <vector>
#include <zlib.h>

namespace kmat {

namespace {

std::string to_lower_ascii(std::string s) {
  for (char& c : s) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return s;
}

bool ends_with_ci(const std::string& path, const std::string& suffix) {
  const std::string p = to_lower_ascii(path);
  const std::string s = to_lower_ascii(suffix);
  if (p.size() < s.size()) {
    return false;
  }
  return p.compare(p.size() - s.size(), s.size(), s) == 0;
}

void strip_cr(std::string& line) {
  if (!line.empty() && line.back() == '\r') {
    line.pop_back();
  }
}

class LineReader {
 public:
  explicit LineReader(const std::string& path) : path_(path), gzip_(ends_with_ci(path, ".gz")) {
    if (gzip_) {
      gz_ = gzopen(path.c_str(), "rb");
      if (gz_ != nullptr) {
        gzbuffer(gz_, 1 << 20);
      }
    } else {
      plain_.open(path);
    }
  }

  ~LineReader() {
    if (gz_ != nullptr) {
      gzclose(gz_);
      gz_ = nullptr;
    }
  }

  LineReader(const LineReader&) = delete;
  LineReader& operator=(const LineReader&) = delete;

  bool ok_open() const { return gzip_ ? (gz_ != nullptr) : static_cast<bool>(plain_); }

  Error getline(std::string& line) {
    line.clear();
    if (gzip_) {
      for (;;) {
        const int ch = gzgetc(gz_);
        if (ch < 0) {
          int err = Z_OK;
          const char* msg = gzerror(gz_, &err);
          if (gzeof(gz_) || err == Z_OK || err == Z_STREAM_END) {
            eof_ = true;
            strip_cr(line);
            return Error::success();
          }
          return Error::io_error(std::string("gzip read failed: ") + path_ + " (" +
                                 (msg != nullptr ? msg : "unknown") + ")");
        }
        if (ch == '\n') {
          strip_cr(line);
          return Error::success();
        }
        line.push_back(static_cast<char>(ch));
      }
    }

    if (!std::getline(plain_, line)) {
      eof_ = true;
      line.clear();
      if (plain_.bad()) {
        return Error::io_error("failed reading file: " + path_);
      }
      return Error::success();
    }
    strip_cr(line);
    return Error::success();
  }

  bool eof() const { return eof_; }

 private:
  std::string path_;
  bool gzip_{false};
  bool eof_{false};
  gzFile gz_{nullptr};
  std::ifstream plain_;
};

Error for_each_fastq_sequence(LineReader& in, const std::string& path,
                              const std::function<Error(const std::string&)>& fn) {
  std::string header;
  std::string seq;
  std::string plus;
  std::string qual;
  std::size_t nrec = 0;

  for (;;) {
    if (auto err = in.getline(header); !err.ok()) {
      return err;
    }
    while (!in.eof() && header.empty()) {
      if (auto err = in.getline(header); !err.ok()) {
        return err;
      }
    }
    if (header.empty()) {
      break;
    }
    if (header.front() != '@') {
      return Error::invalid_argument("FASTQ expected '@' header in " + path);
    }

    if (auto err = in.getline(seq); !err.ok()) {
      return err;
    }
    if (in.eof() && seq.empty()) {
      return Error::invalid_argument("truncated FASTQ record in " + path);
    }
    if (auto err = in.getline(plus); !err.ok()) {
      return err;
    }
    if (auto err = in.getline(qual); !err.ok()) {
      return err;
    }
    if (seq.size() != qual.size()) {
      return Error::invalid_argument("FASTQ sequence/quality length mismatch in " + path);
    }
    if (auto err = normalize_dna(seq); !err.ok()) {
      return err;
    }
    if (auto err = fn(seq); !err.ok()) {
      return err;
    }
    ++nrec;
  }

  if (nrec == 0) {
    return Error::invalid_argument("no sequences found in FASTQ: " + path);
  }
  return Error::success();
}

Error for_each_fasta_sequence(LineReader& in, const std::string& path,
                              const std::function<Error(const std::string&)>& fn) {
  std::string line;
  std::string seq;
  bool have = false;
  std::size_t nrec = 0;

  auto flush = [&]() -> Error {
    if (!have) {
      return Error::success();
    }
    if (auto err = normalize_dna(seq); !err.ok()) {
      return err;
    }
    if (auto err = fn(seq); !err.ok()) {
      return err;
    }
    ++nrec;
    seq.clear();
    have = false;
    return Error::success();
  };

  for (;;) {
    if (auto err = in.getline(line); !err.ok()) {
      return err;
    }
    if (in.eof() && line.empty()) {
      break;
    }
    if (line.empty()) {
      continue;
    }
    if (line.front() == '>') {
      if (auto err = flush(); !err.ok()) {
        return err;
      }
      have = true;
      continue;
    }
    if (!have) {
      return Error::invalid_argument("FASTA sequence before header in " + path);
    }
    seq += line;
  }
  if (auto err = flush(); !err.ok()) {
    return err;
  }
  if (nrec == 0) {
    return Error::invalid_argument("no sequences found in FASTA: " + path);
  }
  return Error::success();
}

}  // namespace

bool path_looks_gzip(const std::string& path) {
  return ends_with_ci(path, ".gz");
}

bool path_looks_fastq(const std::string& path) {
  return ends_with_ci(path, ".fastq") || ends_with_ci(path, ".fq") ||
         ends_with_ci(path, ".fastq.gz") || ends_with_ci(path, ".fq.gz");
}

bool path_looks_fasta(const std::string& path) {
  return ends_with_ci(path, ".fasta") || ends_with_ci(path, ".fa") || ends_with_ci(path, ".fna") ||
         ends_with_ci(path, ".fasta.gz") || ends_with_ci(path, ".fa.gz") ||
         ends_with_ci(path, ".fna.gz");
}

Error for_each_sequence(const std::string& path,
                        const std::function<Error(const std::string& sequence)>& fn) {
  LineReader in(path);
  if (!in.ok_open()) {
    return Error::io_error("failed to open file: " + path);
  }
  if (path_looks_fastq(path)) {
    return for_each_fastq_sequence(in, path, fn);
  }
  return for_each_fasta_sequence(in, path, fn);
}

Error read_fastq_sequences(const std::string& path, std::vector<FastqRecord>& records) {
  records.clear();
  // Re-stream with ids for API compatibility.
  LineReader in(path);
  if (!in.ok_open()) {
    return Error::io_error("failed to open FASTQ: " + path);
  }

  std::string header;
  std::string seq;
  std::string plus;
  std::string qual;

  for (;;) {
    if (auto err = in.getline(header); !err.ok()) {
      return err;
    }
    while (!in.eof() && header.empty()) {
      if (auto err = in.getline(header); !err.ok()) {
        return err;
      }
    }
    if (header.empty()) {
      break;
    }
    if (header.front() != '@') {
      return Error::invalid_argument("FASTQ expected '@' header in " + path);
    }
    if (auto err = in.getline(seq); !err.ok()) {
      return err;
    }
    if (auto err = in.getline(plus); !err.ok()) {
      return err;
    }
    if (auto err = in.getline(qual); !err.ok()) {
      return err;
    }
    if (seq.size() != qual.size()) {
      return Error::invalid_argument("FASTQ sequence/quality length mismatch in " + path);
    }
    FastqRecord rec;
    rec.id = header.size() > 1 ? header.substr(1) : std::string{};
    rec.sequence = seq;
    rec.quality = std::move(qual);
    if (auto err = normalize_dna(rec.sequence); !err.ok()) {
      return err;
    }
    records.push_back(std::move(rec));
  }

  if (records.empty()) {
    return Error::invalid_argument("no sequences found in FASTQ: " + path);
  }
  return Error::success();
}

Error read_sequence_file(const std::string& path, std::vector<FastaRecord>& records) {
  records.clear();
  return for_each_sequence(path, [&](const std::string& sequence) -> Error {
    FastaRecord fr;
    fr.sequence = sequence;
    records.push_back(std::move(fr));
    return Error::success();
  });
}

}  // namespace kmat
