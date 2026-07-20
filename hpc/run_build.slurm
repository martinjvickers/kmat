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
#
# Memory: pass ~90% of the Slurm allocation as --memory-gb so the streaming
# builder can size partition shards. Override with MEMORY_GB=... if needed.

set -euo pipefail

KMAT_IMG="${KMAT_IMG:-$HOME/bin/kmat.img}"
LIST="${LIST:-paths.txt}"
K=31

mkdir -p logs

# Prefer node-local scratch for partition spills when available.
if [ -n "${SLURM_TMPDIR:-}" ]; then
  export TMPDIR="${SLURM_TMPDIR}"
elif [ -d /tmp ]; then
  export TMPDIR="${TMPDIR:-/tmp}"
fi

# Derive --memory-gb from Slurm --mem (MiB) unless MEMORY_GB is set.
if [ -z "${MEMORY_GB:-}" ]; then
  MEM_MB="${SLURM_MEM_PER_NODE:-0}"
  if [ "${MEM_MB}" -gt 0 ] 2>/dev/null; then
    # 90% of allocation, at least 1 GiB
    MEMORY_GB=$(( MEM_MB * 90 / 100 / 1024 ))
    if [ "${MEMORY_GB}" -lt 1 ]; then
      MEMORY_GB=1
    fi
  else
    MEMORY_GB=64
  fi
fi

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
  echo "$(pwd)/${KSET}" >> kset_list.txt
done < "$LIST"

if [ "$missing" -ne 0 ]; then
  echo "ERROR: $missing kset(s) missing — not building"
  exit 1
fi

echo "Building panel.kmat from $(wc -l < kset_list.txt) accessions (memory_gb=${MEMORY_GB})"
singularity exec "$KMAT_IMG" kmat --profile hpc --threads "${SLURM_CPUS_PER_TASK}" \
  build -k "$(pwd)/kset_list.txt" -s "$K" -o "$(pwd)/panel.kmat" \
  --memory-gb "${MEMORY_GB}" --tmpdir "${TMPDIR}"

echo "Done: $(pwd)/panel.kmat"
