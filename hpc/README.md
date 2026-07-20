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

KMC reads `.fq.gz` natively. The image must ship **upstream** KMC under `/usr/local/bin` (not the Debian `/usr/bin/kmc` package, which often fails with “Some error while reading gzip file”).

```bash
sudo singularity build ~/bin/kmat.img singularity/kmat.def
# verify:
singularity exec ~/bin/kmat.img which kmc    # expect /usr/local/bin/kmc
```
