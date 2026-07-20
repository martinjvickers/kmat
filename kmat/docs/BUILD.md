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

Production build:

1. **Scatter** (parallel over accessions): hash-partition each `.kset` by k-mer code into \(T\) sorted buckets under `--tmpdir`
2. **Merge+dedup** (parallel over partitions): multiway-merge the \(N\) buckets for each code partition, build presence words, dedupe patterns → `.pat` / `.map`
3. **Assemble** v2 `.kmat` (pattern concat + wave-merge of sorted maps)

Peak RAM is sized from `--memory-gb`, not from panel unique-kmer count. `--threads` drives both scatter and merge+dedup concurrency.

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
| `--batch-rows` | Requested rows per I/O flush (default **100000**). Scatter batch may be lowered so \(T\) bucket buffers fit. |
| `--tmpdir` | Spill directory for scatter buckets / shards (deleted after success). Prefer `$SLURM_TMPDIR`. |
| `--threads` | Parallel scatter over accessions **and** parallel per-partition merge+dedup. |

Slurm helper [hpc/run_build.slurm](../../hpc/run_build.slurm) sets `set -euo pipefail`, passes `--memory-gb` ≈ 90% of `--mem`, raises `ulimit -n`, and uses `$SLURM_TMPDIR` when present.

## Progress logging

Long phases emit heartbeats about every **30s** (override with `KMAT_BUILD_LOG_EVERY_SEC`):

```text
[kmat:info] scatter: start accessions=440 partitions=64 threads=16
[kmat:info] scatter: 120/440 (27.3%) elapsed=300s rate=0.40/s codes=... eta=800s
[kmat:info] merge_dedup: 8/64 (12.5%) elapsed=45s rate=0.18/s eta=311s
[kmat:info] assemble: 2100000000/8400000000 (25.0%) elapsed=90s rate=23.3e6/s eta=270s
```

Each phase also logs a final `… done` line with totals and elapsed time.

## Scalability model (N=10³–10⁴, large genomes)

| Resource | Bound | Behaviour |
|---|---|---|
| Open files during merge+dedup | `RLIMIT_NOFILE` (soft raised to hard) | Per partition: if live bucket count ≤ budget → one multiway merge; else **hierarchical** reduce. Concurrent workers capped by FD budget. |
| Partition count \(T\) | `max(threads, memory sizing)`, cap **65536** | Code-hash partitions. \(\sum K_i\) upper-bounds unique \(U\) (larger \(T\) = safer RAM). |
| Scatter FDs | open-append-close per bucket flush | Large \(T\) does not hold \(T\) FDs open. |
| Map assemble | FD budget | Patterns concatenated one file at a time; maps wave-merged. |
| Pattern ids | `uint32` | Fails clearly if a shard or total exceeds \(2^{32}-1\). |
| Pattern store size | — | Same presence pattern may appear in multiple code partitions (duplicate rows). Presence per k-mer is correct; GWAS may re-test identical patterns. |

Debug / tests:

- `KMAT_BUILD_MAX_OPEN=<n>` — tiny merge FD budget (hierarchical path)
- `KMAT_BUILD_LOG_EVERY_SEC=<n>` — heartbeat interval

## Memory model

`--memory-gb` budgets per-partition pattern tables (~`memory/threads` each) and scatter buffers. Raising it lowers \(T\) (fewer, larger partitions). Lowering it raises \(T\).

What it does **not** mean: keeping the full unique-kmer table or all accession lists resident.

## Recommended Slurm

| Panel | `--mem` | `--cpus-per-task` | Notes |
|---|---|---|---|
| Toy / medium fixtures | 8–32G | 4–8 | Fine for CI-sized `.kset` |
| ~440 accession teff-scale | 64–256G | 8–16 | Start at 16 threads; prefer node-local scratch |
| 1k–10k accessions / wheat-like \(U\) | 256–512G+ | 16–32 | `$SLURM_TMPDIR`; hierarchical merge if soft `nofile` stays ~1024 |

Do **not** jump to `--threads 64` on NFS scratch without checking I/O: scatter is write-heavy.

## Lessons learned (do not regress)

1. **Stock / apt KMC gzip is unreliable.** The Singularity image builds **patched KMC 3.2.4** (system zlib + `gzread` in `binary_reader`). See [singularity/patches/](../../singularity/patches/).
2. **Never gunzip whole FASTQs to feed KMC**; never load a whole FASTQ into RAM for count.
3. **`.kset` is the handoff**, not KMC DBs: presence-only, sorted `uint64`, KMC-free build.
4. **`gzip -t` (full file), not `zcat \| head`**, is the integrity check for `.fq.gz`.
5. **Singularity `%post` scratch must use `/opt/tmp`**, never host `/tmp` (shared, sticky leftovers).
6. **Build must not hold \(\sum\) accession k-mer lists** or a monolithic `map<code, bitvector>` for production N.
7. **Legacy win was bounded I/O batches (100k rows)** + external parallelism — port that idea; do not treat create-blank-then-fill as the only design.
8. **Slurm scripts need `set -e`** so OOM/killed builds are not reported as success after a trailing `echo Done`.
9. **File-descriptor limits matter, but do not hard-cap partitions for FDs.** Scatter/merge must not hold \(T\) or \(N\cdot T\) FDs. Hierarchical merge when \(N\) exceeds the FD budget.
10. **Keep `singularity exec … kmat … build` as one shell command.** Commenting out only the `singularity exec \` line leaves a bare `build` → `command not found`.
11. **Never size \(T\) from \(\sum K_i\) while also opening \(T\) files.** \(\sum K_i\) as a \(U\) upper bound is fine for *memory* sizing only.
12. **A global single-threaded multiway merge of all `.kset` files will sit at 1× CPU for hours** with no progress logs — partition by code and heartbeat regularly.

## Sequence-file build

`kmat build` still accepts FASTA/FASTQ/`.gz` for small tests. That path loads k-mers in RAM and logs a warning. Production panels must use `count` → `.kset` → memory-bounded build.
