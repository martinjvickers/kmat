# Benchmarks

Harness: `kmat_bench` (built when `KMAT_BUILD_BENCHES=ON`, default).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target kmat_bench

./build/benches/kmat_bench --testdata testdata --panel medium --profile laptop --json /tmp/laptop.json
./build/benches/kmat_bench --testdata testdata --panel medium --profile hpc --json /tmp/hpc.json
```

Panels: `tiny` (k=3, N=6) or `medium` (`panel_k31_n72`, k=31, N=72).

## Runtime profiles

| Profile | Default threads | I/O buffer |
|---|---|---|
| `laptop` | `max(1, hw/2)` | 1 MiB |
| `hpc` | `hardware_concurrency` | 8 MiB |

CLI: `kmat --profile laptop|hpc --threads N <subcommand> …`

## What is timed

1. **build** — sequence → v2 matrix (parallel accession ingest)
2. **pop** — PCA covariates
3. **gwas** — pattern-parallel association (phenotypes residualized once)

## Reference numbers

Recorded on **Darwin arm64** (Apple Silicon), Release build, medium panel. These fixtures are tiny; treat values as smoke baselines for regression detection, not wheat-scale capacity claims.

| Profile | Threads | build_ms | pop_ms | gwas_ms |
|---|---:|---:|---:|---:|
| laptop | 4 | ~5 | ~4 | ~3 |
| hpc | 8 | ~4 | ~5 | ~4 |

Re-run after perf changes; if a step grows by a large factor on the same machine/panel, investigate before merging.

CTest smoke: `kmat_bench_laptop_tiny`, `kmat_bench_hpc_tiny` (must complete successfully).

## Optimisations included in Phase 4

- Laptop / HPC runtime profiles (`--profile`, `--threads`)
- Parallel accession ingest in `build`
- Parallel pattern scoring in `gwas` (deterministic vs 1-thread)
- Rolling k-mer encode (`for_each_encoded_kmer`)
- O(K) pattern→k-mer index for GWAS hit expansion
- Deferred `pa_bits` unless `-d/--display`
- Larger sequential write buffers on HPC profile
