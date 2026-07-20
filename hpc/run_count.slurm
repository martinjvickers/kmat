#!/bin/bash
#SBATCH --job-name=kmat-count
#SBATCH --output=logs/count_%A_%a.out
#SBATCH --error=logs/count_%A_%a.err
#SBATCH --array=1-440%50
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#
# Count k-mers for each accession in paths.txt (one array task per line).
#
# Usage (from the directory that contains paths.txt):
#   mkdir -p logs ksets
#   N=$(grep -cve '^\s*#' -e '^\s*$' paths.txt)   # or: wc -l < paths.txt
#   sbatch --array=1-${N}%50 run_count.slurm

KMAT_IMG="${KMAT_IMG:-$HOME/bin/kmat.img}"
LIST="${LIST:-paths.txt}"
K=31
CI=2

mkdir -p logs ksets

FASTQ=$(grep -vE '^\s*(#|$)' "$LIST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ -z "$FASTQ" ]; then
  echo "No path for array task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

STEM=$(basename "$FASTQ")
STEM=${STEM%.fq.gz}
STEM=${STEM%.fastq.gz}
STEM=${STEM%.fq}
STEM=${STEM%.fastq}

OUT="ksets/${STEM}.kset"
echo "task=${SLURM_ARRAY_TASK_ID}  $FASTQ  ->  $OUT"

singularity exec "$KMAT_IMG" kmat count \
  -i "$FASTQ" -s "$K" --ci "$CI" -o "$OUT"
