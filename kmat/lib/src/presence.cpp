#include "kmat/presence.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>

namespace kmat {

namespace {

bool ends_with_ci(const std::string& path, const std::string& suffix) {
  if (path.size() < suffix.size()) {
    return false;
  }
  for (std::size_t i = 0; i < suffix.size(); ++i) {
    const char a = static_cast<char>(std::tolower(static_cast<unsigned char>(
        path[path.size() - suffix.size() + i])));
    const char b = static_cast<char>(std::tolower(static_cast<unsigned char>(suffix[i])));
    if (a != b) {
      return false;
    }
  }
  return true;
}

}  // namespace

bool path_looks_presence_set(const std::string& path) {
  return ends_with_ci(path, ".kset");
}

Error write_presence_set(const std::string& path, const PresenceSet& set) {
  if (set.header.num_kmers != set.kmers.size()) {
    return Error::invalid_argument("presence set header count mismatch");
  }
  for (std::size_t i = 1; i < set.kmers.size(); ++i) {
    if (set.kmers[i] <= set.kmers[i - 1]) {
      return Error::invalid_argument("presence set k-mers must be strictly sorted");
    }
  }

  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to open presence set for writing: " + path);
  }
  PresenceHeader header = set.header;
  header.magic[0] = 'K';
  header.magic[1] = 'S';
  header.magic[2] = 'E';
  header.magic[3] = 'T';
  header.version = 1;
  header.num_kmers = set.kmers.size();

  out.write(reinterpret_cast<const char*>(&header), sizeof(PresenceHeader));
  if (!set.kmers.empty()) {
    out.write(reinterpret_cast<const char*>(set.kmers.data()),
              static_cast<std::streamsize>(set.kmers.size() * sizeof(std::uint64_t)));
  }
  if (!out) {
    return Error::io_error("failed while writing presence set: " + path);
  }
  return Error::success();
}

Error read_presence_set(const std::string& path, PresenceSet& set) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return Error::io_error("failed to open presence set: " + path);
  }
  set = PresenceSet{};
  in.read(reinterpret_cast<char*>(&set.header), sizeof(PresenceHeader));
  if (!in) {
    return Error::io_error("failed to read presence set header");
  }
  if (set.header.magic[0] != 'K' || set.header.magic[1] != 'S' || set.header.magic[2] != 'E' ||
      set.header.magic[3] != 'T') {
    return Error::invalid_argument("invalid presence set magic (expected KSET)");
  }
  if (set.header.version != 1) {
    return Error::invalid_argument("unsupported presence set version");
  }
  set.kmers.resize(static_cast<std::size_t>(set.header.num_kmers));
  if (!set.kmers.empty()) {
    in.read(reinterpret_cast<char*>(set.kmers.data()),
            static_cast<std::streamsize>(set.kmers.size() * sizeof(std::uint64_t)));
    if (!in) {
      return Error::io_error("unexpected EOF reading presence set k-mers");
    }
  }
  for (std::size_t i = 1; i < set.kmers.size(); ++i) {
    if (set.kmers[i] <= set.kmers[i - 1]) {
      return Error::invalid_argument("presence set k-mers are not sorted: " + path);
    }
  }
  return Error::success();
}

Error PresenceSetCursor::open(const std::string& path) {
  close();
  path_ = path;
  in_.open(path, std::ios::binary);
  if (!in_) {
    return Error::io_error("failed to open presence set: " + path);
  }
  in_.read(reinterpret_cast<char*>(&header_), sizeof(PresenceHeader));
  if (!in_) {
    return Error::io_error("failed to read presence set header: " + path);
  }
  if (header_.magic[0] != 'K' || header_.magic[1] != 'S' || header_.magic[2] != 'E' ||
      header_.magic[3] != 'T') {
    return Error::invalid_argument("invalid presence set magic (expected KSET): " + path);
  }
  if (header_.version != 1) {
    return Error::invalid_argument("unsupported presence set version: " + path);
  }
  remaining_ = header_.num_kmers;
  ok_ = true;
  has_value_ = false;
  return advance();
}

void PresenceSetCursor::close() {
  in_.close();
  ok_ = false;
  has_value_ = false;
  remaining_ = 0;
  path_.clear();
}

Error PresenceSetCursor::advance() {
  has_value_ = false;
  if (!ok_) {
    return Error::io_error("presence set cursor is not open");
  }
  if (remaining_ == 0) {
    return Error::success();
  }
  const std::uint64_t prev = value_;
  in_.read(reinterpret_cast<char*>(&value_), sizeof(std::uint64_t));
  if (!in_) {
    ok_ = false;
    return Error::io_error("unexpected EOF reading presence set k-mers: " + path_);
  }
  if (remaining_ < header_.num_kmers && value_ <= prev) {
    ok_ = false;
    return Error::invalid_argument("presence set k-mers are not sorted: " + path_);
  }
  --remaining_;
  has_value_ = true;
  return Error::success();
}

}  // namespace kmat
