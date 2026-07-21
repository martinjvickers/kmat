#!/bin/bash
#SBATCH --job-name=kmat-fill
#SBATCH --output=logs/fill_%A_%a.out
#SBATCH --error=logs/fill_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=0-0
#
# Fill one stripe's columns from .kset files (sequential -n loop per task).
# Requires stripes/list_XX.txt with ≤64 .kset paths for that stripe.
# Set --array to cover all stripes.

set -euo pipefail

KMAT_IMG="${KMAT_IMG:-$HOME/bin/kmat.img}"
PREFIX="${PREFIX:-panel}"
IDX="${SLURM_ARRAY_TASK_ID:-0}"
STRIPE=$(printf "%s.%02d.bin" "$PREFIX" "$IDX")
LIST=$(printf "stripes/list_%02d.txt" "$IDX")

mkdir -p logs

if [ ! -f "stripes/${STRIPE}" ]; then
  echo "ERROR: missing stripes/${STRIPE}"
  exit 1
fi
if [ ! -f "$LIST" ]; then
  echo "ERROR: missing $LIST"
  exit 1
fi

WORKDIR="${PWD}"
FINAL_STRIPE="${WORKDIR}/stripes/${STRIPE}"
STAGE_STRIPE="$FINAL_STRIPE"
STAGED=0

# Stage to node-local scratch only when it is a different path (avoids mv same-file noise).
if [ -n "${SLURM_LOCAL_SCRATCH:-}" ]; then
  STAGE_DIR="${SLURM_LOCAL_SCRATCH}/kmat_fill_${SLURM_JOB_ID:-$$}_${IDX}"
  mkdir -p "$STAGE_DIR"
  cp "$FINAL_STRIPE" "${STAGE_DIR}/${STRIPE}"
  STAGE_STRIPE="${STAGE_DIR}/${STRIPE}"
  STAGED=1
fi

n=0
while IFS= read -r KSET; do
  case "$KSET" in ''|\#*) continue ;; esac
  echo "Filling n=$n from $KSET"
  singularity exec "$KMAT_IMG" kmat --profile hpc fill \
    -o "$STAGE_STRIPE" --kset "$KSET" -n "$n" --batch-rows 100000
  n=$((n + 1))
done < "${WORKDIR}/${LIST}"

if [ "$STAGED" -eq 1 ]; then
  mv "$STAGE_STRIPE" "$FINAL_STRIPE"
  rmdir "$(dirname "$STAGE_STRIPE")" 2>/dev/null || true
fi
echo "Done: stripes/${STRIPE} columns=$n"
