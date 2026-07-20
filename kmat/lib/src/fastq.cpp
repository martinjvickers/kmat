#include "kmat/fastq.hpp"

#include "kmat/fasta.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
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

Error read_all_bytes_plain(const std::string& path, std::string& out) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return Error::io_error("failed to open file: " + path);
  }
  out.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  return Error::success();
}

Error read_all_bytes_gzip(const std::string& path, std::string& out) {
  gzFile gz = gzopen(path.c_str(), "rb");
  if (gz == nullptr) {
    return Error::io_error("failed to open gzip file: " + path);
  }

  out.clear();
  char buf[1 << 16];
  while (true) {
    const int n = gzread(gz, buf, sizeof(buf));
    if (n < 0) {
      gzclose(gz);
      return Error::io_error("gzip read failed: " + path);
    }
    if (n == 0) {
      break;
    }
    out.append(buf, static_cast<std::size_t>(n));
  }
  gzclose(gz);
  return Error::success();
}

Error read_all_bytes(const std::string& path, std::string& out) {
  if (path_looks_gzip(path)) {
    return read_all_bytes_gzip(path, out);
  }
  return read_all_bytes_plain(path, out);
}

std::vector<std::string> split_lines(const std::string& text) {
  std::vector<std::string> lines;
  std::string cur;
  for (char c : text) {
    if (c == '\n') {
      if (!cur.empty() && cur.back() == '\r') {
        cur.pop_back();
      }
      lines.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  if (!cur.empty()) {
    if (cur.back() == '\r') {
      cur.pop_back();
    }
    lines.push_back(cur);
  }
  return lines;
}

Error parse_fastq_text(const std::string& text, const std::string& path,
                       std::vector<FastqRecord>& records) {
  const std::vector<std::string> lines = split_lines(text);
  records.clear();

  for (std::size_t i = 0; i < lines.size();) {
    while (i < lines.size() && lines[i].empty()) {
      ++i;
    }
    if (i >= lines.size()) {
      break;
    }
    if (lines[i].empty() || lines[i].front() != '@') {
      return Error::invalid_argument("FASTQ expected '@' header in " + path);
    }
    if (i + 3 >= lines.size()) {
      return Error::invalid_argument("truncated FASTQ record in " + path);
    }

    FastqRecord rec;
    rec.id = lines[i].substr(1);
    rec.sequence = lines[i + 1];
    // lines[i + 2] is '+' separator (optional id)
    rec.quality = lines[i + 3];
    if (rec.sequence.size() != rec.quality.size()) {
      return Error::invalid_argument("FASTQ sequence/quality length mismatch in " + path);
    }
    if (auto err = normalize_dna(rec.sequence); !err.ok()) {
      return err;
    }
    records.push_back(std::move(rec));
    i += 4;
  }

  if (records.empty()) {
    return Error::invalid_argument("no sequences found in FASTQ: " + path);
  }
  return Error::success();
}

Error parse_fasta_text(const std::string& text, const std::string& path,
                       std::vector<FastaRecord>& records) {
  const std::vector<std::string> lines = split_lines(text);
  records.clear();
  FastaRecord current;
  bool have_record = false;

  for (const std::string& line : lines) {
    if (line.empty()) {
      continue;
    }
    if (line.front() == '>') {
      if (have_record) {
        if (auto err = normalize_dna(current.sequence); !err.ok()) {
          return err;
        }
        records.push_back(std::move(current));
        current = FastaRecord{};
      }
      current.id = line.substr(1);
      have_record = true;
      continue;
    }
    if (!have_record) {
      return Error::invalid_argument("FASTA sequence before header in " + path);
    }
    current.sequence += line;
  }

  if (have_record) {
    if (auto err = normalize_dna(current.sequence); !err.ok()) {
      return err;
    }
    records.push_back(std::move(current));
  }

  if (records.empty()) {
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

Error read_fastq_sequences(const std::string& path, std::vector<FastqRecord>& records) {
  std::string text;
  if (auto err = read_all_bytes(path, text); !err.ok()) {
    return err;
  }
  return parse_fastq_text(text, path, records);
}

Error read_sequence_file(const std::string& path, std::vector<FastaRecord>& records) {
  if (path_looks_fastq(path)) {
    std::vector<FastqRecord> fq;
    if (auto err = read_fastq_sequences(path, fq); !err.ok()) {
      return err;
    }
    records.clear();
    records.reserve(fq.size());
    for (FastqRecord& r : fq) {
      FastaRecord fr;
      fr.id = std::move(r.id);
      fr.sequence = std::move(r.sequence);
      records.push_back(std::move(fr));
    }
    return Error::success();
  }

  // FASTA (plain or gzip) or unknown → try FASTA text parse after optional gunzip.
  if (path_looks_gzip(path) || path_looks_fasta(path)) {
    std::string text;
    if (auto err = read_all_bytes(path, text); !err.ok()) {
      return err;
    }
    return parse_fasta_text(text, path, records);
  }

  // Default: existing plain FASTA reader.
  return read_fasta_sequences(path, records);
}

}  // namespace kmat
