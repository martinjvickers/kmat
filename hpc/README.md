# Running kmat on HPC

This folder has a Singularity image recipe and Slurm helpers for a **count-then-build** panel (e.g. ~100 `.fastq.gz` accessions).

## 1. Build the Singularity image

```bash
sudo singularity build kmat.img singularity/kmat.def
```

Needs network (apt + clone `https://github.com/martinjvickers/kmat.git` + CMake deps).


## 2. Prepare an accession list

One path per line (absolute paths that will be visible inside the container via `--bind`):

```text
/project/wheat/fastq/acc_001.fastq.gz
/project/wheat/fastq/acc_002.fastq.gz
...
```

See [`accession_list.example.txt`](accession_list.example.txt). Blank lines and `#` comments are ignored. **List order becomes matrix column order** (same IDs GWAS/pop will use after stripping extensions).

## 3. Submit count + build

```bash
cd /scratch/$USER/kmat_runs   # or wherever you keep work
cp /path/to/kmat.sif .
cp /path/to/my_100.txt ./accession_list.txt
mkdir -p logs

export KMAT_SIF=$PWD/kmat.sif
export ACCESSION_LIST=$PWD/accession_list.txt
export OUTDIR=$PWD/out_panel100
export BIND=/project:/project,/scratch:/scratch   # site-specific
export KMER_SIZE=31
export CI=2                    # drop low-abundance k-mers (KMC-style -ci)

# from the cloned repo:
/path/to/kmat_repo/hpc/submit_panel.sh
```

What this does:

1. **Array job** — one task per line in `ACCESSION_LIST`:
   `kmat count -i <fastq.gz> -s 31 --ci 2 -o out/ksets/<id>.kset`
2. **Build job** (after array succeeds) — writes `out/panel.kmat` from those `.kset` files in list order.

Or submit the two steps yourself:

```bash
N=$(grep -vE '^\s*(#|$)' accession_list.txt | wc -l | tr -d ' ')
COUNT=$(sbatch --parsable --array=1-${N}%50 \
  --export=ALL,KMAT_SIF=$PWD/kmat.sif,ACCESSION_LIST=$PWD/accession_list.txt,OUTDIR=$PWD/out,CI=2,BIND=/project:/project \
  /path/to/kmat_repo/hpc/slurm_count_array.sh)

sbatch --dependency=afterok:${COUNT} \
  --export=ALL,KMAT_SIF=$PWD/kmat.sif,ACCESSION_LIST=$PWD/accession_list.txt,OUTDIR=$PWD/out,BIND=/project:/project \
  /path/to/kmat_repo/hpc/slurm_build.sh
```

## 4. After the matrix exists

Use the **kset list** (same column order) for validate/pop/gwas:

```bash
# built automatically as ${OUTDIR}/kset_list.txt
singularity exec --bind "$BIND" "$KMAT_SIF" \
  kmat validate -i "$OUTDIR/panel.kmat" -k "$OUTDIR/kset_list.txt"

singularity exec --bind "$BIND" "$KMAT_SIF" \
  kmat --profile hpc pop -i "$OUTDIR/panel.kmat" -k "$OUTDIR/kset_list.txt" -o "$OUTDIR/pop.tsv"

# phenotypes.tsv: accession IDs must match stems (acc_001 from acc_001.fastq.gz)
singularity exec --bind "$BIND" "$KMAT_SIF" \
  kmat --profile hpc gwas -i "$OUTDIR/panel.kmat" -k "$OUTDIR/kset_list.txt" \
    -p phenotypes.tsv --pop "$OUTDIR/pop.tsv" -s 31 -N 1000
```

## Memory warning (important)

Current `kmat count` loads each FASTQ into memory before counting. For large accessions, raise `#SBATCH --mem` in `slurm_count_array.sh` (start at 32–64G and watch `seff`). Streaming ingest is a follow-up improvement if you hit OOMs.

## Tuning

| Variable | Meaning |
|---|---|
| `CI` | Min k-mer count (`--ci`) |
| `KMER_SIZE` | `-s` (default 31) |
| `ARRAY_THROTTLE` | Max concurrent array tasks (`%50` in submit script) |
| `KMAT_PROFILE` | `hpc` (default) or `laptop` |
| `BIND` | Singularity/Apptainer bind mounts for data paths |

Edit `#SBATCH` headers in the scripts for your partition/account/QoS.
