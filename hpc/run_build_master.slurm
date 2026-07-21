#!/bin/bash
#SBATCH --job-name=kmat-master
#SBATCH --output=logs/master_%j.out
#SBATCH --error=logs/master_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#
# Tree-merge ksets/ → panel.kuniv (few intermediate files, not N×T).
# Parallelism = min(threads, ceil(N/group_size)) independent group merges.
# Usage after count array:
#   sbatch --dependency=afterok:<COUNT_JOBID> run_build_master.slurm

set -euo pipefail

KMAT_IMG="${KMAT_IMG:-$HOME/bin/kmat.img}"
LIST="${LIST:-paths.txt}"
K=31
GROUP="${GROUP_SIZE:-32}"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

mkdir -p logs
if [ -n "${SLURM_TMPDIR:-}" ]; then
  export TMPDIR="${SLURM_TMPDIR}"
fi
mkdir -p "${TMPDIR:-/tmp}"

: > kset_list.txt
missing=0
while IFS= read -r FASTQ; do
  case "$FASTQ" in ''|\#*) continue ;; esac
  STEM=$(basename "$FASTQ")
  STEM=${STEM%.fq.gz}; STEM=${STEM%.fastq.gz}; STEM=${STEM%.fq}; STEM=${STEM%.fastq}
  KSET="ksets/${STEM}.kset"
  if [ ! -f "$KSET" ]; then
    echo "MISSING $KSET"
    missing=$((missing + 1))
    continue
  fi
  echo "$(pwd)/${KSET}" >> kset_list.txt
done < "$LIST"

if [ "$missing" -ne 0 ]; then
  echo "ERROR: $missing kset(s) missing"
  exit 1
fi

N=$(wc -l < kset_list.txt)
echo "Building panel.kuniv from ${N} accessions (threads=${THREADS} group_size=${GROUP})"
singularity exec "$KMAT_IMG" kmat --profile hpc --threads "$THREADS" build-master \
  -k "$(pwd)/kset_list.txt" -s "$K" -o "$(pwd)/panel.kuniv" \
  --tmpdir "${TMPDIR:-/tmp}" --group-size "$GROUP"

echo "Done: $(pwd)/panel.kuniv"
