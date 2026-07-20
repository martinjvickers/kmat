#!/usr/bin/env python3
"""Generate the medium k=31 / N=72 FASTQ.gz toy panel (GitHub-sized).

Re-run from this directory:
  python3 generate.py

Outputs are committed so CI does not need to regenerate.
"""

from __future__ import annotations

import gzip
import pathlib
import random

N_ACC = 72
K = 31
SEED = 42
OUT = pathlib.Path(__file__).resolve().parent

# Shared background long enough for many 31-mers; deterministic.
_rng = random.Random(SEED)
BASES = "ACGT"


def rand_dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def write_fastq_gz(path: pathlib.Path, read_id: str, seq: str) -> None:
    qual = "I" * len(seq)
    payload = f"@{read_id}\n{seq}\n+\n{qual}\n".encode("ascii")
    with gzip.open(path, "wb", compresslevel=9) as fh:
        fh.write(payload)


def main() -> None:
    acc_dir = OUT / "accessions"
    acc_dir.mkdir(parents=True, exist_ok=True)

    shared = rand_dna(360, _rng)
    # Planted marker present only in accessions 0..35 (first half) for GWAS signal.
    marker = "ACGTACGTACGTACGTACGTACGTACGTA"  # 29 bp; with flanks becomes unique 31-mers
    marker_block = "AA" + marker + "TT"

    list_lines: list[str] = []
    pheno_lines = ["accession\ttrait1"]

    for i in range(N_ACC):
        name = f"acc_{i:03d}"
        path = acc_dir / f"{name}.fastq.gz"
        unique = rand_dna(45, random.Random(SEED + 1000 + i))
        seq = shared + unique
        if i < 36:
            seq = seq + marker_block
        write_fastq_gz(path, name, seq)
        list_lines.append(f"accessions/{name}.fastq.gz")
        # High trait for marker carriers, low otherwise, plus tiny noise.
        trait = (10.0 if i < 36 else 0.0) + (i % 5) * 0.01
        pheno_lines.append(f"{name}\t{trait:.4f}")

    (OUT / "accession_list.txt").write_text("\n".join(list_lines) + "\n", encoding="utf-8")
    (OUT / "phenotypes.tsv").write_text("\n".join(pheno_lines) + "\n", encoding="utf-8")

    # Gene spans start of shared region so gene search finds shared 31-mers.
    gene_seq = shared[:80]
    (OUT / "gene.fasta").write_text(f">shared_gene\n{gene_seq}\n", encoding="utf-8")

    print(f"Wrote {N_ACC} fastq.gz accessions under {acc_dir}")
    print(f"Shared length={len(shared)} gene_len={len(gene_seq)} k={K}")


if __name__ == "__main__":
    main()
