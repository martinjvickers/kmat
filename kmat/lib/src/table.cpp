#include "kmat/table.hpp"

#include <fstream>
#include <sstream>

namespace kmat {

namespace {

std::vector<std::string> split_tab(const std::string& line) {
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (start <= line.size()) {
    const std::size_t pos = line.find('\t', start);
    if (pos == std::string::npos) {
      fields.push_back(line.substr(start));
      break;
    }
    fields.push_back(line.substr(start, pos - start));
    start = pos + 1;
  }
  return fields;
}

}  // namespace

Error read_numeric_table(const std::string& path, std::vector<std::string>& header,
                         std::vector<TableRow>& rows) {
  std::ifstream in(path);
  if (!in) {
    return Error::io_error("failed to open table: " + path);
  }

  header.clear();
  rows.clear();
  std::string line;
  bool got_header = false;
  std::size_t expected_cols = 0;

  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line.empty() || line.front() == '#') {
      continue;
    }

    const auto fields = split_tab(line);
    if (fields.empty()) {
      continue;
    }

    if (!got_header) {
      header = fields;
      if (header.size() < 2) {
        return Error::invalid_argument("table must have accession column plus data columns");
      }
      expected_cols = header.size();
      got_header = true;
      continue;
    }

    if (fields.size() != expected_cols) {
      return Error::invalid_argument("inconsistent column count in " + path);
    }

    TableRow row;
    row.key = fields[0];
    row.values.reserve(fields.size() - 1);
    for (std::size_t i = 1; i < fields.size(); ++i) {
      try {
        row.values.push_back(std::stod(fields[i]));
      } catch (const std::exception&) {
        return Error::invalid_argument("non-numeric value in " + path + " row " + row.key);
      }
    }
    rows.push_back(std::move(row));
  }

  if (!got_header || rows.empty()) {
    return Error::invalid_argument("table is empty: " + path);
  }
  return Error::success();
}

Error write_pop_tsv(const std::string& path, const std::vector<std::string>& accessions,
                    const std::vector<std::vector<double>>& pcs) {
  if (accessions.size() != pcs.size()) {
    return Error::invalid_argument("accession count does not match PC row count");
  }
  if (accessions.empty()) {
    return Error::invalid_argument("no accessions to write");
  }

  const std::size_t npc = pcs.front().size();
  for (const auto& row : pcs) {
    if (row.size() != npc) {
      return Error::invalid_argument("inconsistent PC column count");
    }
  }

  std::ofstream out(path, std::ios::trunc);
  if (!out) {
    return Error::io_error("failed to open pop TSV for writing: " + path);
  }

  out << "accession";
  for (std::size_t j = 0; j < npc; ++j) {
    out << '\t' << "PC" << (j + 1);
  }
  out << '\n';

  for (std::size_t i = 0; i < accessions.size(); ++i) {
    out << accessions[i];
    for (std::size_t j = 0; j < npc; ++j) {
      out << '\t' << pcs[i][j];
    }
    out << '\n';
  }

  if (!out) {
    return Error::io_error("failed while writing pop TSV: " + path);
  }
  return Error::success();
}

}  // namespace kmat
