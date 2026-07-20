# `kmat build` — memory-bounded matrix construction

## Pipeline

```text
kmat count  →  *.kset   (per accession, Slurm array)
kmat build  →  panel.kmat  (one job, streaming union)
```

`.kset` files are sorted unique `uint64` presence codes (see [FORMAT.md](FORMAT.md)). Build never re-reads FASTQ and never opens KMC databases.

## Why streaming

The first in-RAM builder loaded every `.kset` into `vector`s and inserted all unique k-mers into a `std::map<code, bitvector>`. For N≈440 that peaks at roughly:

\[
8\cdot\sum_i K_i \;+\; O(U\cdot S)
\]

bytes (plus allocator / tree overhead), which OOMs even on large HPC nodes.

Production build **multiway-merges** sorted `.kset` streams, **hash-partitions** by presence-pattern fingerprint into spill shards under `--tmpdir`, **deduplicates patterns per shard** (parallel), then assembles a v2 `.kmat`. Peak RAM is sized from `--memory-gb`, not from panel unique-kmer count.

## CLI

```bash
kmat --profile hpc --threads 16 build \
  -k kset_list.txt -s 31 -o panel.kmat \
  --memory-gb 64 \
  --batch-rows 100000 \
  --tmpdir "${TMPDIR:-/tmp}"
```

| Flag | Meaning |
|---|---|
| `--memory-gb` | Working-set budget (GiB). Drives partition count. `0` → profile default (laptop **8**, hpc **64**). |
| `--batch-rows` | Rows per shard I/O flush (default **100000**, same quantum as legacy fill). |
| `--tmpdir` | Spill directory for partition shards (deleted after success). |
| `--threads` | Parallel shard pattern-dedup (and global runtime workers). |

Slurm helper [hpc/run_build.slurm](../../hpc/run_build.slurm) sets `set -euo pipefail`, passes `--memory-gb` ≈ 90% of `--mem`, and uses `$SLURM_TMPDIR` when present.

## Memory model

`--memory-gb` budgets **one shard worker’s** pattern table + row buffers (roughly `memory / threads` per parallel worker). Partition count \(P\) is chosen from \(\sum K_i\) (upper bound on unique \(U\)) so each shard’s expected working set fits. Raising `--memory-gb` lowers \(P\) (fewer, larger shards). Lowering it raises \(P\) (more spill, less RAM).

What it does **not** mean: keeping the full unique-kmer table or all accession lists resident.

## Recommended Slurm

| Panel | `--mem` | `--cpus-per-task` | Notes |
|---|---|---|---|
| Toy / medium fixtures | 8–32G | 4–8 | Fine for CI-sized `.kset` |
| ~440 accession teff-scale | 64–256G | 8–16 | Start at 64G; raise if a shard still OOMs |
| Wheat-like unique counts | 256–512G | 16–32 | More RAM → fewer partitions → less I/O |

Do **not** set `--threads` near 64 unless `--mem` is large: more threads means more concurrent shard tables.

## Lessons learned (do not regress)

1. **Stock / apt KMC gzip is unreliable.** The Singularity image builds **patched KMC 3.2.4** (system zlib + `gzread` in `binary_reader`). See [singularity/patches/](../../singularity/patches/).
2. **Never gunzip whole FASTQs to feed KMC**; never load a whole FASTQ into RAM for count.
3. **`.kset` is the handoff**, not KMC DBs: presence-only, sorted `uint64`, KMC-free build.
4. **`gzip -t` (full file), not `zcat \| head`**, is the integrity check for `.fq.gz`.
5. **Singularity `%post` scratch must use `/opt/tmp`**, never host `/tmp` (shared, sticky leftovers).
6. **Build must not hold \(\sum\) accession k-mer lists** or a monolithic `map<code, bitvector>` for production N.
7. **Legacy win was bounded I/O batches (100k rows)** + external parallelism — port that idea; do not treat create-blank-then-fill as the only design.
8. **Slurm scripts need `set -e`** so OOM/killed builds are not reported as success after a trailing `echo Done`.

## Sequence-file build

`kmat build` still accepts FASTA/FASTQ/`.gz` for small tests. That path loads k-mers in RAM and logs a warning. Production panels must use `count` → `.kset` → memory-bounded build.
