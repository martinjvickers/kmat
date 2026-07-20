#pragma once

#include <string>

namespace kmat {

/// Lightweight status for library operations (no exceptions required).
enum class ErrorCode {
  Ok = 0,
  InvalidArgument,
  IoError,
  InvalidKmer,
};

struct Error {
  ErrorCode code{ErrorCode::Ok};
  std::string message;

  [[nodiscard]] bool ok() const { return code == ErrorCode::Ok; }
  explicit operator bool() const { return !ok(); }

  static Error success() { return {}; }

  static Error invalid_argument(std::string msg) {
    return {ErrorCode::InvalidArgument, std::move(msg)};
  }

  static Error io_error(std::string msg) {
    return {ErrorCode::IoError, std::move(msg)};
  }

  static Error invalid_kmer(std::string msg) {
    return {ErrorCode::InvalidKmer, std::move(msg)};
  }
};

}  // namespace kmat
