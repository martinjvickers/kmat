#pragma once

#include "kmat/error.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace kmat {

struct TableRow {
  std::string key;
  std::vector<double> values;
};

/// Read a tab-delimited table with accession in column 0 and numeric columns after.
Error read_numeric_table(const std::string& path, std::vector<std::string>& header,
                         std::vector<TableRow>& rows);

Error write_pop_tsv(const std::string& path, const std::vector<std::string>& accessions,
                    const std::vector<std::vector<double>>& pcs);

}  // namespace kmat
