#include "kmat/validate.hpp"

#include "kmat/io.hpp"
#include "kmat/matrix.hpp"

namespace kmat {

Error validate_panel(const ValidateOptions& opts, ValidateReport& report) {
  report = ValidateReport{};
  report.ok = true;
  report.summary = "validation passed";

  std::vector<std::string> accession_paths;
  if (auto err = read_list_file(opts.accession_list_path, accession_paths); !err.ok()) {
    report.ok = false;
    report.summary = "failed to read accession list";
    report.issues.push_back(err.message);
    return Error::success();
  }
  if (accession_paths.empty()) {
    report.ok = false;
    report.summary = "accession list is empty";
    report.issues.emplace_back("no accessions");
    return Error::success();
  }

  PaMatrix matrix;
  if (!opts.matrix_list_path.empty()) {
    if (auto err = load_matrix_from_list(opts.matrix_list_path, accession_paths.size(), matrix);
        !err.ok()) {
      report.ok = false;
      report.summary = "failed to load matrix from list";
      report.issues.push_back(err.message);
      return Error::success();
    }
  } else if (!opts.matrix_path.empty()) {
    if (auto err = read_matrix(opts.matrix_path, matrix); !err.ok()) {
      report.ok = false;
      report.summary = "failed to read matrix";
      report.issues.push_back(err.message);
      return Error::success();
    }
  } else {
    report.ok = false;
    report.summary = "matrix path or matrix list required";
    report.issues.emplace_back("no matrix input");
    return Error::success();
  }

  if (matrix.header.num_accessions != accession_paths.size()) {
    report.ok = false;
    report.summary = "accession count mismatch";
    report.issues.push_back("header num_accessions=" + std::to_string(matrix.header.num_accessions) +
                            " list=" + std::to_string(accession_paths.size()));
  }

  const std::size_t expected_stripes = stripe_count_for_accessions(accession_paths.size());
  if (matrix.header.num_stripes != expected_stripes) {
    report.ok = false;
    report.summary = "stripe count mismatch";
    report.issues.push_back("header num_stripes=" + std::to_string(matrix.header.num_stripes) +
                            " expected=" + std::to_string(expected_stripes));
  }

  if (matrix.kmers.empty()) {
    report.ok = false;
    report.summary = "matrix has no k-mers";
    report.issues.emplace_back("num_kmers=0");
  } else if (matrix.header.num_rows != matrix.kmers.size()) {
    report.ok = false;
    report.summary = "k-mer count mismatch";
    report.issues.push_back("header num_rows=" + std::to_string(matrix.header.num_rows) +
                            " actual=" + std::to_string(matrix.kmers.size()));
  }

  if (matrix.header.reserved != matrix.patterns.size()) {
    report.ok = false;
    report.summary = "pattern count mismatch";
    report.issues.push_back("header reserved=" + std::to_string(matrix.header.reserved) +
                            " actual=" + std::to_string(matrix.patterns.size()));
  }

  for (std::size_t i = 1; i < matrix.kmers.size(); ++i) {
    if (matrix.kmers[i].kmer_code <= matrix.kmers[i - 1].kmer_code) {
      report.ok = false;
      report.summary = "k-mer map not strictly sorted";
      report.issues.emplace_back("sort order broken at entry " + std::to_string(i));
      break;
    }
  }

  for (std::size_t i = 0; i < matrix.kmers.size(); ++i) {
    if (matrix.kmers[i].pattern_id >= matrix.patterns.size()) {
      report.ok = false;
      report.summary = "k-mer map pattern_id out of range";
      report.issues.push_back("entry " + std::to_string(i));
      break;
    }
  }

  for (std::size_t i = 0; i < matrix.patterns.size(); ++i) {
    if (matrix.patterns[i].size() != matrix.header.num_stripes) {
      report.ok = false;
      report.summary = "inconsistent pattern stripe width";
      report.issues.push_back("pattern " + std::to_string(i));
      break;
    }
  }

  if (report.ok) {
    report.summary = "matrix OK: " + std::to_string(matrix.kmers.size()) + " k-mers, " +
                     std::to_string(matrix.patterns.size()) + " patterns, " +
                     std::to_string(accession_paths.size()) + " accessions, k=" +
                     std::to_string(matrix.header.kmer_size) + ", v" +
                     std::to_string(matrix.header.version);
  }

  return Error::success();
}

}  // namespace kmat
