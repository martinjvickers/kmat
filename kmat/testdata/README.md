# Synthetic panels for kmat Phase 1

## Tiny panel (default)

| Path | Role |
|---|---|
| `accessions/*.fasta` | Per-accession sequences (k=3 smoke tests) |
| `accession_list.txt` | Accession FASTA paths |
| `phenotypes.tsv` | Phenotype table |
| `gene.fasta` | Short gene for `kmat gene` |

## Medium panel

See [`panel_k31_n72/`](panel_k31_n72/): **k=31**, **72 accessions**, **`.fastq.gz`**, ~300 KB — exercises multi-stripe matrices and gzip FASTQ ingest. Sized for GitHub.

## Matrix format (Phase 1)

Binary `.kmat` file:

- 40-byte header: magic `KMAT`, version, k, num accessions, num stripes, num rows
- Each row: `uint64 kmer_code` + `num_stripes` × `uint64` presence words (LSB = accession 0 within the stripe)

`num_stripes = ceil(N/64)`. For the medium panel, N=72 → 2 stripes.

## Example pipelines

Tiny (from repo `kmat/`):

```bash
TD=testdata
OUT=/tmp/kmat_panel
mkdir -p "$OUT"

./build/cli/kmat build -k "$TD/accession_list.txt" -s 3 -o "$OUT/panel.kmat"
./build/cli/kmat validate -i "$OUT/panel.kmat" -k "$TD/accession_list.txt"
./build/cli/kmat pop -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" -o "$OUT/pop.tsv"
./build/cli/kmat gwas -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" \
  -p "$TD/phenotypes.tsv" --pop "$OUT/pop.tsv" -s 3 -a > "$OUT/gwas.tsv"
./build/cli/kmat gene -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" \
  -g "$TD/gene.fasta" -s 3 > "$OUT/gene.tsv"
```

Medium:

```bash
TD=testdata/panel_k31_n72
OUT=/tmp/kmat_k31_n72
mkdir -p "$OUT"

./build/cli/kmat build -k "$TD/accession_list.txt" -s 31 -o "$OUT/panel.kmat"
./build/cli/kmat validate -i "$OUT/panel.kmat" -k "$TD/accession_list.txt"
# … pop / gwas / gene with -s 31 (see panel_k31_n72/README.md)
```
