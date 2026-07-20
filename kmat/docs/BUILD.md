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
| `--batch-rows` | Requested rows per shard I/O flush (default **100000**). May be lowered automatically so \(P\) writer buffers fit in ~¼ of the budget. |
| `--tmpdir` | Spill directory for partition / merge shards (deleted after success). |
| `--threads` | Parallel shard pattern-dedup (and global runtime workers). |

Slurm helper [hpc/run_build.slurm](../../hpc/run_build.slurm) sets `set -euo pipefail`, passes `--memory-gb` ≈ 90% of `--mem`, and uses `$SLURM_TMPDIR` when present.

## Scalability model (N=10³–10⁴, large genomes)

| Resource | Bound | Behaviour |
|---|---|---|
| Open files during merge | `RLIMIT_NOFILE` (soft raised to hard at start) | If \(N\) ≤ budget → one multiway merge. If \(N\) > budget → **hierarchical** group merge into sorted row spills, then reduce. |
| Partition count \(P\) | Memory + absolute cap **65536** | Sized so expected rows/shard fit in `memory/threads`. \(\sum K_i\) is only an **upper bound** on unique \(U\) (makes \(P\) larger / safer). **Not** capped at 128. |
| Shard FDs while writing | 1 at a time | Writers **open-append-close** on flush — large \(P\) does not mean large concurrent FDs. |
| Map assemble | FD budget | Pattern shards concatenated one-at-a-time; map shards multiway-merged in **waves** if \(P\) is large. |
| Pattern ids | `uint32` | Build fails clearly if a shard or the total exceeds \(2^{32}-1\). |

Do **not** treat teff (N≈440) as the design point: wheat-scale panels need hierarchical merge and memory-driven \(P\), not hard FD×partition coupling.

Debug / tests: `KMAT_BUILD_MAX_OPEN=<n>` forces a tiny merge FD budget (exercises hierarchical path).

## Memory model

`--memory-gb` budgets per-shard pattern tables (~`memory/threads` each) and aggregate writer buffers (~¼ budget). Raising it lowers \(P\) (fewer, larger shards). Lowering it raises \(P\) (more spill, less RAM per shard).

What it does **not** mean: keeping the full unique-kmer table or all accession lists resident.

## Recommended Slurm

| Panel | `--mem` | `--cpus-per-task` | Notes |
|---|---|---|---|
| Toy / medium fixtures | 8–32G | 4–8 | Fine for CI-sized `.kset` |
| ~440 accession teff-scale | 64–256G | 8–16 | Start at 64G; raise if a shard still OOMs |
| 1k–10k accessions / wheat-like \(U\) | 256–512G+ | 16–32 | Prefer node-local `$SLURM_TMPDIR`; hierarchical merge if soft `nofile` stays ~1024 |

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
9. **File-descriptor limits matter, but do not hard-cap partitions for FDs.** Merge may keep one FD per live input stream; shard writers must not hold \(P\) FDs. When \(N\) exceeds the FD budget, merge hierarchically. A mistaken `P≤128` hard cap makes large-\(U\) shards OOM.
10. **Keep `singularity exec … kmat … build` as one shell command.** Commenting out only the `singularity exec \` line leaves a bare `build` → `command not found`.
11. **Never size \(P\) from \(\sum K_i\) while also opening \(P\) files.** \(\sum K_i\) as a \(U\) upper bound is fine for *memory* sizing only.

## Sequence-file build

`kmat build` still accepts FASTA/FASTQ/`.gz` for small tests. That path loads k-mers in RAM and logs a warning. Production panels must use `count` → `.kset` → memory-bounded build.
