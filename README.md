# kmat (refactor workspace)

Greenfield **kmat** toolkit plus roadmap docs. Legacy `kmer_search/` / `matrix_presetup/` are gitignored reference trees.

| Path | Role |
|---|---|
| [`kmat/`](kmat/) | Product source (CLI + lib + tests) |
| [`hpc/`](hpc/) | Singularity + Slurm helpers for real panels |
| [`singularity/kmat.def`](singularity/kmat.def) | Container recipe |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | Legacy pipeline reference |
| [`REFACTOR.md`](REFACTOR.md) | Target design + phase tracker |

## Quick start (laptop)

```bash
cd kmat
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## HPC

```bash
sudo singularity build ~/bin/kmat.img singularity/kmat.def   # includes KMC + kmat
# then see hpc/README.md — count uses KMC (--engine kmc)
```
