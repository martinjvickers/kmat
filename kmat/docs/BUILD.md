# `kmat build` — multi-node, few-file matrix construction

## Pipeline

```text
kmat count          →  *.kset          (Slurm array, per accession)
kmat build-master   →  panel.kuniv     (tree-merge sorted .kset → unique codes)
kmat create-stripe  →  panel.XX.bin    (Slurm array over stripes; blank dense)
kmat fill           →  panel.XX.bin    (set bits from .kset; 100k-row batches)
kmat compress       →  panel.kmat      (v2 pattern-compressed; global dedup)
```

Small panels: `kmat build -k kset_list.txt …` runs these stages in-process (still few files under `--tmpdir`).

Production panels: use the HPC scripts under [`hpc/`](../../hpc/) (array jobs + local SSD), same shape as Watkins `matrix_presetup/`.

## Why this shape

Legacy production ([`ARCHITECTURE.md`](../../ARCHITECTURE.md)) won scale with:

- **Few large files** — O(N/64) stripe bins, not hundreds of thousands of spills
- **Column parallelism** — Slurm array over stripes across nodes
- **100k-row batch I/O** — fill never loads U×N
- **Local SSD staging** — rewrite traffic off NFS

An earlier scatter/hash-shard design (`a{i}/p{t}.bin`) fixed single-node OOM but **regressed** those properties (N·T inodes, single-node only, duplicate patterns). It is **not** the production path.

## Durable file budget

| Artifact | Count |
|---|---|
| `.kuniv` master | **1** |
| Dense stripe `.bin` | **ceil(N/64)** |
| Final v2 `.kmat` | **1** |
| Tree-merge temps | O(N/G) then deleted (G≈16–64) |

## CLI

```bash
# Master universe (sorted unique codes; parallel group merges)
kmat --profile hpc --threads 32 build-master -k kset_list.txt -s 31 -o panel.kuniv \
  --tmpdir "$TMPDIR" --group-size 32

# Blank stripe (≤64 accessions in stripe list)
kmat create-stripe -m panel.kuniv -k list_00.txt -s 31 -o panel.00.bin \
  --num-accessions 351 --stripe-index 0 --batch-rows 100000

# Fill one local column
kmat fill -k list_00.txt -n 7 -o panel.00.bin --batch-rows 100000

# Compress filled stripes → v2
kmat compress -m matrix_list.txt -k unified_list.txt -s 31 -o panel.kmat --memory-gb 64
```

| Flag / stage | Meaning |
|---|---|
| `--threads` | **build-master**: parallel independent group merges (and reduce). Effective workers ≈ `min(threads, ceil(N/group_size))`. |
| `--group-size` | Fan-in per merge node (default 32). Smaller → more parallel groups; larger → fewer temps / FDs. |
| `--batch-rows` | Create/fill I/O window (default **100000**, legacy default). |
| `--memory-gb` | Bounds **compress** pattern-dictionary working set (spill to a small number of large shards if needed). Does **not** invent N×T files. |
| `--tmpdir` | Scratch for master tree-merge temps / in-process build. Prefer `$SLURM_TMPDIR`. |

Progress: **build-master** logs every ~15s / 5M unique codes per group (`unique=… rate=…/s`), plus `groups_done=i/G`. Other stages: every **30s** (`KMAT_BUILD_LOG_EVERY_SEC`).

## HPC

| Script | Role |
|---|---|
| [`hpc/run_build_master.slurm`](../../hpc/run_build_master.slurm) | Tree-merge → `panel.kuniv` |
| [`hpc/run_create_stripes.slurm`](../../hpc/run_create_stripes.slurm) | Array over stripes; SSD stage |
| [`hpc/run_fill.slurm`](../../hpc/run_fill.slurm) | Array over stripes; fill local `-n` loop; SSD stage |
| [`hpc/run_compress.slurm`](../../hpc/run_compress.slurm) | Stripes → `panel.kmat` |
| [`hpc/run_build.slurm`](../../hpc/run_build.slurm) | Small-panel convenience: single-node staged build |

## Lessons learned (do not regress)

1. **Stock / apt KMC gzip is unreliable** — Singularity builds patched KMC 3.2.4.
2. **Never gunzip whole FASTQs for KMC**; never load a whole FASTQ into RAM for count.
3. **`.kset` is the handoff**, not KMC DBs, for fill.
4. **`gzip -t`**, not `zcat \| head`, for integrity.
5. **Singularity `%post` scratch: `/opt/tmp`**, not host `/tmp`.
6. **Build must not hold Σ accession k-mer lists** or a monolithic `map<code, bitvector>`.
7. **100k-row batches + stripe arrays + local SSD** are sacred orchestration — not optional.
8. **Slurm scripts need `set -e`**.
9. **Never N×partition tiny spill files** for production.
10. **Keep `singularity exec … kmat …` as one shell command**.
11. **Global pattern dedup at compress** — no intentional cross-shard duplicate patterns.
12. **`--memory-gb` must not invent inode storms**.

## Sequence-file build

`kmat build` still accepts FASTA/FASTQ/`.gz` for small tests (loads k-mers in RAM; warns). Production: `count` → `.kset` → master/create/fill/compress.
