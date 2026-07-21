#!/bin/bash
#SBATCH --job-name=kmat-compress
#SBATCH --output=logs/compress_%j.out
#SBATCH --error=logs/compress_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#
# Compress filled stripes/ → panel.kmat (v2, globally deduped patterns).

set -euo pipefail

KMAT_IMG="${KMAT_IMG:-$HOME/bin/kmat.img}"
K=31
ACC_LIST="${ACC_LIST:-kset_list.txt}"
MATRIX_LIST="${MATRIX_LIST:-stripes/matrix_list.txt}"

mkdir -p logs
if [ -n "${SLURM_TMPDIR:-}" ]; then
  export TMPDIR="${SLURM_TMPDIR}"
fi

MEM_GB="${MEMORY_GB:-}"
if [ -z "$MEM_GB" ]; then
  MEM_MB="${SLURM_MEM_PER_NODE:-0}"
  if [ "${MEM_MB}" -gt 0 ] 2>/dev/null; then
    MEM_GB=$(( MEM_MB * 90 / 100 / 1024 ))
    [ "$MEM_GB" -lt 1 ] && MEM_GB=1
  else
    MEM_GB=64
  fi
fi

singularity exec "$KMAT_IMG" kmat --profile hpc compress \
  -m "$(pwd)/${MATRIX_LIST}" -k "$(pwd)/${ACC_LIST}" -s "$K" \
  -o "$(pwd)/panel.kmat" --memory-gb "$MEM_GB" --batch-rows 100000 \
  --tmpdir "${TMPDIR:-/tmp}"

echo "Done: $(pwd)/panel.kmat"
