# Simple Slurm usage (no binds)

Production `kmat count` uses **KMC** inside the Singularity image (`--engine kmc`).

## Small panels (single-node staged build)

```bash
cd /path/to/testing_kmat
mkdir -p logs ksets

N=$(grep -vE '^\s*(#|$)' paths.txt | wc -l)
COUNT=$(sbatch --parsable --array=1-${N}%50 run_count.slurm)
sbatch --dependency=afterok:${COUNT} run_build.slurm
```

`run_build.slurm` runs master → create/fill stripes → v2 compress in one job (few large files under `$SLURM_TMPDIR`).

## Production panels (multi-node, few-file)

Same shape as Watkins `matrix_presetup/`:

```bash
COUNT=$(sbatch --parsable --array=1-${N}%50 run_count.slurm)
MASTER=$(sbatch --parsable --dependency=afterok:${COUNT} run_build_master.slurm)
# Set --array=0-$((STRIPES-1)) and NACC=<accession count>
CREATE=$(sbatch --parsable --dependency=afterok:${MASTER} --array=0-$((STRIPES-1)) \
  --export=ALL,NACC=${NACC} run_create_stripes.slurm)
FILL=$(sbatch --parsable --dependency=afterok:${CREATE} --array=0-$((STRIPES-1)) run_fill.slurm)
sbatch --dependency=afterok:${FILL} run_compress.slurm
```

| Script | Role |
|---|---|
| `run_count.slurm` | Array → `ksets/*.kset` |
| `run_build_master.slurm` | Tree-merge → `panel.kuniv` |
| `run_create_stripes.slurm` | Array → blank `stripes/panel.XX.bin` (SSD stage when available) |
| `run_fill.slurm` | Array → fill columns from `.kset` (100k-row batches) |
| `run_compress.slurm` | Stripes → `panel.kmat` (v2, global pattern dedup) |
| `run_build.slurm` | Small-panel all-in-one staged build |

Details: [`kmat/docs/BUILD.md`](../kmat/docs/BUILD.md).

Defaults: `$HOME/bin/kmat.img`, `paths.txt`, `-s 31 --ci 2`.

KMC reads `.fq.gz` natively. The image builds a **patched KMC 3.2.4** (system zlib + `gzread`).

```bash
sudo singularity build ~/bin/kmat.img singularity/kmat.def
singularity exec ~/bin/kmat.img which kmc
```
