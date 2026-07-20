#pragma once

#include "kmat/error.hpp"

#include <string>
#include <vector>

namespace kmat {

/// Read a newline-delimited list file, skipping blank lines and lines starting with '#'.
Error read_list_file(const std::string& path, std::vector<std::string>& lines);

Error write_list_file(const std::string& path, const std::vector<std::string>& lines);

/// Make relative entries in a list file absolute using the list file's directory.
Error resolve_list_paths(const std::string& list_path, std::vector<std::string>& paths);

/// Basename without compound extensions (.fastq.gz, .fa.gz, .fasta, …).
std::string accession_id_from_path(const std::string& path);

}  // namespace kmat
