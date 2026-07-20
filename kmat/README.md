# kmat

Greenfield rebuild of the k-mer presence/absence GWAS toolkit. Phase 4 adds laptop/HPC runtime profiles, parallel build/GWAS paths, and a `kmat_bench` harness with recorded baselines.

For **HPC / Singularity / Slurm** (count ~100 FASTQ.gz → matrix), see [`../hpc/README.md`](../hpc/README.md).

## Requirements

- CMake 3.16+
- C++17 compiler (Apple Clang, GCC, or Clang)
- zlib (system / Homebrew)
- Git (FetchContent: CLI11, Catch2, Eigen)
- For production `count --engine kmc`: `kmc` and `kmc_tools` on `PATH` (bundled in Singularity image)

## Build (macOS / Linux)

```bash
cd kmat
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Disable tests / benches:

```bash
cmake -S . -B build -DKMAT_BUILD_TESTS=OFF -DKMAT_BUILD_BENCHES=OFF
cmake --build build
```

## Runtime profiles

```bash
./build/cli/kmat --profile laptop --threads 4 build ...
./build/cli/kmat --profile hpc gwas ...
```

| Profile | Threads (default) | Matrix write buffer |
|---|---|---|
| `laptop` | `hw/2` | 1 MiB |
| `hpc` | all cores | 8 MiB |

See [`docs/BENCH.md`](docs/BENCH.md).

## Pipelines

### Tiny panel (k=3, 6 accessions, FASTA)

From the `kmat/` directory:

```bash
TD=testdata
OUT=/tmp/kmat_panel
mkdir -p "$OUT"

./build/cli/kmat build -k "$TD/accession_list.txt" -s 3 -o "$OUT/panel.kmat"
./build/cli/kmat validate -i "$OUT/panel.kmat" -k "$TD/accession_list.txt"
./build/cli/kmat pop -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" -o "$OUT/pop.tsv"
./build/cli/kmat gwas -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" \
  -p "$TD/phenotypes.tsv" --pop "$OUT/pop.tsv" -s 3 -a
./build/cli/kmat gene -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" \
  -g "$TD/gene.fasta" -s 3
```

### Count then build (HPC-friendly)

Production path uses **KMC** (`--engine kmc`, default). Needs `kmc` + `kmc_tools` on `PATH` (Singularity image includes them). For tiny local tests without KMC, use `--engine builtin`.

```bash
TD=testdata
OUT=/tmp/kmat_count
mkdir -p "$OUT/ksets"

./build/cli/kmat --threads 8 count --engine builtin \
  -i "$TD/accessions/sample_a.fasta" -s 3 --ci 1 -o "$OUT/ksets/sample_a.kset"
# production / Singularity:
# ./build/cli/kmat --profile hpc --threads 8 count --engine kmc --tmpdir "$TMPDIR" \
#   -i acc.fq.gz -s 31 --ci 2 -o acc.kset

./build/cli/kmat build -k "$OUT/kset_list.txt" -s 3 -o "$OUT/panel.kmat"
```

### Medium panel (k=31, 72 accessions, FASTQ.gz)

```bash
TD=testdata/panel_k31_n72
OUT=/tmp/kmat_k31_n72
mkdir -p "$OUT"

./build/cli/kmat --profile hpc build -k "$TD/accession_list.txt" -s 31 -o "$OUT/panel.kmat"
./build/cli/kmat validate -i "$OUT/panel.kmat" -k "$TD/accession_list.txt"
./build/cli/kmat pop -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" -o "$OUT/pop.tsv"
./build/cli/kmat --profile hpc gwas -i "$OUT/panel.kmat" -k "$TD/accession_list.txt" \
  -p "$TD/phenotypes.tsv" --pop "$OUT/pop.tsv" -s 31 -a
```

## Subcommands

| Command | Status |
|---|---|
| `count` | FASTQ/FASTA → `.kset` via **KMC** (`--engine kmc`) or builtin hashmap |
| `import-kmers` | Text k-mer list → `.kset` (no KMC link) |
| `build` | PA matrix from sequences or `.kset` (parallel ingest) |
| `pop` | PCA population-structure TSV |
| `gwas` | Pattern-parallel association |
| `gene` | Gene FASTA k-mer lookup |
| `validate` | Matrix/list consistency checks |
| `fill` | Stub |

## Benchmarks

```bash
./build/benches/kmat_bench --testdata testdata --panel medium --profile laptop
./build/benches/kmat_bench --testdata testdata --panel medium --profile hpc --json out.json
```

## Test

```bash
cd build
ctest --output-on-failure
```

## Docs

- [`docs/FORMAT.md`](docs/FORMAT.md) — `.kmat` / `.kset` layouts
- [`docs/BENCH.md`](docs/BENCH.md) — profiles and baselines
- [`../REFACTOR.md`](../REFACTOR.md) — roadmap
- [`../ARCHITECTURE.md`](../ARCHITECTURE.md) — legacy reference
