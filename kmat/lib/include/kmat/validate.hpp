#pragma once

#include "kmat/error.hpp"

#include <string>
#include <vector>

namespace kmat {

struct ValidateOptions {
  std::string matrix_path;
  std::string matrix_list_path;
  std::string accession_list_path;
};

struct ValidateReport {
  bool ok{false};
  std::string summary;
  std::vector<std::string> issues;
};

Error validate_panel(const ValidateOptions& opts, ValidateReport& report);

}  // namespace kmat
