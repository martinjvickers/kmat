#!/bin/bash
#SBATCH --job-name=kmat-build
#SBATCH --output=logs/build_%j.out
#SBATCH --error=logs/build_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#
# Build panel.kmat from ksets/ produced by run_count.slurm.
# Uses the same order as paths.txt.
#
# Usage (after the count array finishes):
#   sbatch --dependency=afterok:<COUNT_JOBID> run_build.slurm

KMAT_IMG="${KMAT_IMG:-$HOME/bin/kmat.img}"
LIST="${LIST:-paths.txt}"
K=31

mkdir -p logs

# Build a kset list in the same order as paths.txt
: > kset_list.txt
missing=0
while IFS= read -r FASTQ; do
  case "$FASTQ" in ''|\#*) continue ;; esac
  STEM=$(basename "$FASTQ")
  STEM=${STEM%.fq.gz}
  STEM=${STEM%.fastq.gz}
  STEM=${STEM%.fq}
  STEM=${STEM%.fastq}
  KSET="ksets/${STEM}.kset"
  if [ ! -f "$KSET" ]; then
    echo "MISSING $KSET"
    missing=$((missing + 1))
    continue
  fi
  # Prefer absolute paths so singularity cwd does not matter
  echo "$(pwd)/${KSET}" >> kset_list.txt
done < "$LIST"

if [ "$missing" -ne 0 ]; then
  echo "ERROR: $missing kset(s) missing — not building"
  exit 1
fi

echo "Building panel.kmat from $(wc -l < kset_list.txt) accessions"
singularity exec "$KMAT_IMG" kmat --profile hpc --threads "${SLURM_CPUS_PER_TASK}" \
  build -k "$(pwd)/kset_list.txt" -s "$K" -o "$(pwd)/panel.kmat"

echo "Done: $(pwd)/panel.kmat"
