#!/bin/bash
#SBATCH --job-name=kmat-build
#SBATCH --output=logs/build_%j.out
#SBATCH --error=logs/build_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#
# Small-panel convenience: single-node staged build
#   .kset list → master → create/fill stripes → v2 panel.kmat
#
# For production multi-node panels prefer:
#   run_build_master.slurm → run_create_stripes.slurm → run_fill.slurm → run_compress.slurm
#
# Usage (after the count array finishes):
#   sbatch --dependency=afterok:<COUNT_JOBID> run_build.slurm

set -euo pipefail

KMAT_IMG="${KMAT_IMG:-$HOME/bin/kmat.img}"
LIST="${LIST:-paths.txt}"
K=31
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p logs

if [ -n "${SLURM_TMPDIR:-}" ]; then
  export TMPDIR="${SLURM_TMPDIR}"
elif [ -d /tmp ]; then
  export TMPDIR="${TMPDIR:-/tmp}"
fi
mkdir -p "${TMPDIR}"

if [ -z "${MEMORY_GB:-}" ]; then
  MEM_MB="${SLURM_MEM_PER_NODE:-0}"
  if [ "${MEM_MB}" -gt 0 ] 2>/dev/null; then
    MEMORY_GB=$(( MEM_MB * 90 / 100 / 1024 ))
    if [ "${MEMORY_GB}" -lt 1 ]; then
      MEMORY_GB=1
    fi
  else
    MEMORY_GB=64
  fi
fi

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
  echo "$(pwd)/${KSET}" >> kset_list.txt
done < "$LIST"

if [ "$missing" -ne 0 ]; then
  echo "ERROR: $missing kset(s) missing — not building"
  exit 1
fi

echo "Staged build panel.kmat from $(wc -l < kset_list.txt) accessions (memory_gb=${MEMORY_GB} threads=${THREADS})"
singularity exec "$KMAT_IMG" kmat --profile hpc --threads "${THREADS}" \
  build -k "$(pwd)/kset_list.txt" -s "$K" -o "$(pwd)/panel.kmat" \
  --memory-gb "${MEMORY_GB}" --batch-rows 100000 --tmpdir "${TMPDIR}"

echo "Done: $(pwd)/panel.kmat"
