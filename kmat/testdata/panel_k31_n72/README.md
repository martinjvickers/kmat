# Medium panel: k=31, 72 accessions, FASTQ.gz

GitHub-sized fixture for multi-stripe matrices and gzip FASTQ ingest.

| Item | Value |
|---|---|
| Accessions | 72 (`ceil(72/64)=2` stripes) |
| k | 31 |
| Inputs | `accessions/acc_XXX.fastq.gz` |
| Size | ~300 KB total |

## Files

| Path | Role |
|---|---|
| `accessions/*.fastq.gz` | One short read per accession |
| `accession_list.txt` | Relative paths to FASTQ.gz files |
| `phenotypes.tsv` | Synthetic trait correlated with a planted marker (acc 0–35) |
| `gene.fasta` | Slice of the shared background for `kmat gene` |
| `generate.py` | Reproducible generator (fixtures are committed) |

## Regenerate

```bash
python3 generate.py
```

## Example pipeline (from `kmat/`)

```bash
TD=testdata/panel_k31_n72
OUT=/tmp/kmat_k31_n72
mkdir -p "$OUT"

./build/cli/kmat build -k "$TD/accession_list.txt" -s 31 -o "$OUT/panel.kmat"
./build/cli/kmat validate -i "$OUT/panel.kmat" -k "$TD/accession_list.txt"
./build/cli/kmat pop -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" -o "$OUT/pop.tsv"
./build/cli/kmat gwas -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" \
  -p "$TD/phenotypes.tsv" --pop "$OUT/pop.tsv" -s 31 -a > "$OUT/gwas.tsv"
./build/cli/kmat gene -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" \
  -g "$TD/gene.fasta" -s 31 > "$OUT/gene.tsv"
```
