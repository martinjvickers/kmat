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
# anywhere — def is self-contained (clones kmat for patches)
sudo singularity build ~/bin/kmat.img singularity/kmat.def
# or after copying the def to $HOME as def.def:
#   sudo singularity build ~/bin/kmat.img def.def

singularity exec ~/bin/kmat.img which kmc          # /usr/local/bin/kmc
singularity exec ~/bin/kmat.img ldd /usr/local/bin/kmc   # should list libz.so
```
