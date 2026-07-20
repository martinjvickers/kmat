#include "kmat/gene.hpp"

#include "kmat/fasta.hpp"
#include "kmat/io.hpp"
#include "kmat/kmer.hpp"
#include "kmat/matrix.hpp"
#include "kmat/sequence.hpp"

#include <algorithm>
#include <unordered_set>

namespace kmat {

Error run_gene_search(const GeneSearchOptions& opts, GeneSearchResult& result) {
  std::vector<std::string> accession_paths;
  if (auto err = read_list_file(opts.accession_list_path, accession_paths); !err.ok()) {
    return err;
  }

  PaMatrix matrix;
  if (!opts.matrix_list_path.empty()) {
    if (auto err = load_matrix_from_list(opts.matrix_list_path, accession_paths.size(), matrix);
        !err.ok()) {
      return err;
    }
  } else if (!opts.matrix_path.empty()) {
    if (auto err = read_matrix(opts.matrix_path, matrix); !err.ok()) {
      return err;
    }
  } else {
    return Error::invalid_argument("matrix path or matrix list required");
  }

  if (matrix.header.kmer_size != opts.kmer_size) {
    return Error::invalid_argument("k-mer size mismatch between CLI and matrix");
  }

  std::vector<FastaRecord> records;
  if (auto err = read_fasta_sequences(opts.gene_fasta_path, records); !err.ok()) {
    return err;
  }

  std::unordered_set<std::uint64_t> query_kmers;
  for (const FastaRecord& rec : records) {
    if (auto err = gene_kmer_set(rec.sequence, opts.kmer_size, query_kmers); !err.ok()) {
      return err;
    }
  }

  result.hits.clear();
  std::vector<std::uint64_t> sorted_query(query_kmers.begin(), query_kmers.end());
  std::sort(sorted_query.begin(), sorted_query.end());

  for (std::uint64_t code : sorted_query) {
    const auto pid = find_pattern_id(matrix, code);
    if (!pid) {
      continue;
    }
    GeneHit hit;
    hit.kmer_code = code;
    if (auto err = decode_kmer(code, opts.kmer_size, hit.kmer); !err.ok()) {
      return err;
    }
    hit.pa_bits = presence_bitstring(matrix.patterns[*pid], accession_paths.size());
    result.hits.push_back(std::move(hit));
  }

  return Error::success();
}

Error write_gene_hits(std::ostream& out, const GeneSearchResult& result) {
  for (const GeneHit& hit : result.hits) {
    out << hit.kmer << '\t' << hit.pa_bits << '\n';
  }
  return Error::success();
}

}  // namespace kmat
