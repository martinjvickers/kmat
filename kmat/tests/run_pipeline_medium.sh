#!/usr/bin/env bash
set -euo pipefail

KMAT_BIN="${1:?kmat binary path required}"
TD="${2:?panel_k31_n72 directory required}"
OUT="$(mktemp -d)"

cleanup() { rm -rf "$OUT"; }
trap cleanup EXIT

"$KMAT_BIN" build -k "$TD/accession_list.txt" -s 31 -o "$OUT/panel.kmat"
"$KMAT_BIN" validate -i "$OUT/panel.kmat" -k "$TD/accession_list.txt"
"$KMAT_BIN" pop -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" -o "$OUT/pop.tsv"
"$KMAT_BIN" gwas -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" \
  -p "$TD/phenotypes.tsv" --pop "$OUT/pop.tsv" -s 31 -a > "$OUT/gwas.tsv"
test -s "$OUT/gwas.tsv"
"$KMAT_BIN" gene -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" \
  -g "$TD/gene.fasta" -s 31 > "$OUT/gene.tsv"
test -s "$OUT/gene.tsv"

echo "medium pipeline OK"
