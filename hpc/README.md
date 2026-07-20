# Simple Slurm usage (no binds, no repo scripts required at runtime)

Copy these two files next to your `paths.txt`:

- `run_count.slurm` — one array task per accession → `ksets/*.kset`
- `run_build.slurm` — union into `panel.kmat`

```bash
cd /path/to/testing_kmat          # has paths.txt
mkdir -p logs ksets

# copy the .slurm files here, then:
N=$(grep -vE '^\s*(#|$)' paths.txt | wc -l)
sbatch --array=1-${N}%50 run_count.slurm

# when that job finishes (note the job id from sbatch):
sbatch --dependency=afterok:JOBID run_build.slurm
```

Defaults: image `$HOME/bin/kmat.img`, list `paths.txt`, `-s 31 --ci 2`.
Override if needed: `sbatch --export=ALL,KMAT_IMG=...,LIST=... ...`
