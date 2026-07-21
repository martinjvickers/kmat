#!/bin/bash
#SBATCH --job-name=kmat-create
#SBATCH --output=logs/create_%A_%a.out
#SBATCH --error=logs/create_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-0
#
# Create blank dense stripes from panel.kuniv (one array task per stripe).
# Set --array=0-$((NUM_STRIPES-1)) to match ceil(N/64).
#
# Prefers local scratch when available (legacy Watkins pattern).

set -euo pipefail

KMAT_IMG="${KMAT_IMG:-$HOME/bin/kmat.img}"
K=31
UNIV="${UNIV:-panel.kuniv}"
NACC="${NACC:-}"
PREFIX="${PREFIX:-panel}"

mkdir -p logs stripes

if [ -z "$NACC" ]; then
  echo "ERROR: set NACC to global accession count"
  exit 1
fi

IDX="${SLURM_ARRAY_TASK_ID:-0}"
OUT_NAME=$(printf "%s.%02d.bin" "$PREFIX" "$IDX")

WORKDIR="${PWD}"
mkdir -p "${WORKDIR}/stripes"
STAGE_DIR="${WORKDIR}"
UNIV_PATH="${WORKDIR}/${UNIV}"

if [ -n "${SLURM_LOCAL_SCRATCH:-}" ]; then
  STAGE_DIR="${SLURM_LOCAL_SCRATCH}/kmat_create_${SLURM_JOB_ID:-$$}_${IDX}"
  mkdir -p "$STAGE_DIR"
  cp "${WORKDIR}/${UNIV}" "${STAGE_DIR}/"
  UNIV_PATH="${STAGE_DIR}/$(basename "$UNIV")"
elif [ -n "${SLURM_TMPDIR:-}" ]; then
  STAGE_DIR="${SLURM_TMPDIR}/kmat_create_${SLURM_JOB_ID:-$$}_${IDX}"
  mkdir -p "$STAGE_DIR"
  cp "${WORKDIR}/${UNIV}" "${STAGE_DIR}/"
  UNIV_PATH="${STAGE_DIR}/$(basename "$UNIV")"
fi

singularity exec "$KMAT_IMG" kmat --profile hpc create-stripe \
  -m "$UNIV_PATH" -s "$K" -o "${STAGE_DIR}/${OUT_NAME}" \
  --num-accessions "$NACC" --stripe-index "$IDX" --batch-rows 100000

if [ "${STAGE_DIR}/${OUT_NAME}" != "${WORKDIR}/stripes/${OUT_NAME}" ]; then
  mv "${STAGE_DIR}/${OUT_NAME}" "${WORKDIR}/stripes/${OUT_NAME}"
  rm -f "${STAGE_DIR}/$(basename "$UNIV")" 2>/dev/null || true
  rmdir "$STAGE_DIR" 2>/dev/null || true
fi
echo "Done: ${WORKDIR}/stripes/${OUT_NAME}"
