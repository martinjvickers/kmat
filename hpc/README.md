# Simple Slurm usage (no binds)

Production `kmat count` uses **KMC** inside the Singularity image (`--engine kmc`).
Copy these next to `paths.txt`:

- `run_count.slurm` — array → `ksets/*.kset` (KMC, multithreaded)
- `run_build.slurm` — union → `panel.kmat`

```bash
cd /path/to/testing_kmat
mkdir -p logs ksets

N=$(grep -vE '^\s*(#|$)' paths.txt | wc -l)
COUNT=$(sbatch --parsable --array=1-${N}%50 run_count.slurm)
sbatch --dependency=afterok:${COUNT} run_build.slurm
```

Defaults: `$HOME/bin/kmat.img`, `paths.txt`, `-s 31 --ci 2`, `--cpus-per-task=8` for KMC `-t`.

KMC reads `.fq.gz` natively. The image builds a **patched KMC 3.2.4** (system zlib + `gzread`) — stock release/apt KMC often fails with “Some error while reading gzip file”.

```bash
# from the repo root (so %files can find singularity/patches/)
sudo singularity build ~/bin/kmat.img singularity/kmat.def
singularity exec ~/bin/kmat.img which kmc          # /usr/local/bin/kmc
singularity exec ~/bin/kmat.img ldd $(which kmc)   # should list libz.so
```
