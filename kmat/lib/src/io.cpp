#include "kmat/io.hpp"

#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace kmat {

namespace {

bool is_ignorable_line(std::string_view line) {
  if (line.empty()) {
    return true;
  }
  return line.front() == '#';
}

}  // namespace

Error read_list_file(const std::string& path, std::vector<std::string>& lines) {
  std::ifstream in(path);
  if (!in) {
    return Error::io_error("failed to open list file for reading: " + path);
  }

  lines.clear();
  std::string raw;
  while (std::getline(in, raw)) {
    if (!raw.empty() && raw.back() == '\r') {
      raw.pop_back();
    }
    if (is_ignorable_line(raw)) {
      continue;
    }
    lines.push_back(raw);
  }

  if (!in.eof()) {
    return Error::io_error("failed while reading list file: " + path);
  }
  return Error::success();
}

Error write_list_file(const std::string& path, const std::vector<std::string>& lines) {
  std::ofstream out(path, std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to open list file for writing: " + path);
  }

  for (const auto& line : lines) {
    out << line << '\n';
  }
  if (!out) {
    return Error::io_error("failed while writing list file: " + path);
  }
  return Error::success();
}

Error resolve_list_paths(const std::string& list_path, std::vector<std::string>& paths) {
  const std::filesystem::path base = std::filesystem::path(list_path).parent_path();
  for (std::string& path : paths) {
    const std::filesystem::path p(path);
    if (p.is_absolute()) {
      continue;
    }
    path = (base / p).lexically_normal().string();
  }
  return Error::success();
}

std::string accession_id_from_path(const std::string& path) {
  const std::size_t slash = path.find_last_of("/\\");
  std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);

  auto ends_with_ci = [](const std::string& s, const std::string& suffix) {
    if (s.size() < suffix.size()) {
      return false;
    }
    for (std::size_t i = 0; i < suffix.size(); ++i) {
      const char a = static_cast<char>(std::tolower(static_cast<unsigned char>(s[s.size() - suffix.size() + i])));
      const char b = static_cast<char>(std::tolower(static_cast<unsigned char>(suffix[i])));
      if (a != b) {
        return false;
      }
    }
    return true;
  };

  static const char* kExts[] = {".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz", ".fna.gz",
                                ".fastq",    ".fq",    ".fasta",    ".fa",    ".fna",
                                ".kset",     ".gz"};
  for (const char* ext : kExts) {
    if (ends_with_ci(base, ext)) {
      base.resize(base.size() - std::char_traits<char>::length(ext));
      break;
    }
  }
  return base;
}

}  // namespace kmat
