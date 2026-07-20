#!/bin/bash
#SBATCH --job-name=kmat-count
#SBATCH --output=logs/kmat_count_%A_%a.out
#SBATCH --error=logs/kmat_count_%A_%a.err
#SBATCH --array=1-100%50
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#
# Per-accession k-mer count → .kset
#
# Required environment (pass with --export=ALL,VAR=... or a small wrapper):
#   KMAT_SIF          path to kmat.sif
#   ACCESSION_LIST    text file: one FASTQ/FASTA (.gz ok) path per line
#   OUTDIR            output directory (ksets/ written here)
# Optional:
#   KMER_SIZE         default 31
#   CI                min count (--ci), default 2
#   BIND              singularity bind specs, e.g. /project:/project,/scratch:/scratch
#   KMAT_PROFILE      laptop|hpc (default hpc)
#
# Submit (example):
#   mkdir -p logs
#   N=$(grep -vE '^\s*(#|$)' "$ACCESSION_LIST" | wc -l | tr -d ' ')
#   sbatch --array=1-${N}%50 \
#     --export=ALL,KMAT_SIF=$PWD/kmat.sif,ACCESSION_LIST=$PWD/my_100.txt,OUTDIR=$PWD/out,CI=2,BIND=/data:/data \
#     hpc/slurm_count_array.sh

set -euo pipefail

: "${KMAT_SIF:?set KMAT_SIF to the kmat singularity image}"
: "${ACCESSION_LIST:?set ACCESSION_LIST to a file of FASTQ paths}"
: "${OUTDIR:?set OUTDIR to the run output directory}"

KMER_SIZE="${KMER_SIZE:-31}"
CI="${CI:-2}"
KMAT_PROFILE="${KMAT_PROFILE:-hpc}"
BIND="${BIND:-}"

mkdir -p "${OUTDIR}/ksets" logs

# Resolve the Nth non-comment, non-blank line (1-based SLURM_ARRAY_TASK_ID).
mapfile -t LINES < <(grep -vE '^\s*(#|$)' "${ACCESSION_LIST}")
N="${#LINES[@]}"
if [[ "${N}" -eq 0 ]]; then
  echo "ERROR: no paths in ${ACCESSION_LIST}" >&2
  exit 1
fi
if [[ "${SLURM_ARRAY_TASK_ID}" -lt 1 || "${SLURM_ARRAY_TASK_ID}" -gt "${N}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range 1..${N}" >&2
  exit 1
fi

FASTQ="${LINES[$((SLURM_ARRAY_TASK_ID - 1))]}"
BASE="$(basename "${FASTQ}")"
# Strip common sequence extensions for the .kset stem.
STEM="${BASE}"
for ext in .fastq.gz .fq.gz .fasta.gz .fa.gz .fna.gz .fastq .fq .fasta .fa .fna .gz; do
  case "${STEM}" in
    *"${ext}") STEM="${STEM%${ext}}"; break ;;
  esac
done
KSET="${OUTDIR}/ksets/${STEM}.kset"

echo "[$(date -Is)] task=${SLURM_ARRAY_TASK_ID}/${N}"
echo "  input : ${FASTQ}"
echo "  output: ${KSET}"
echo "  -s ${KMER_SIZE} --ci ${CI}"

SING_BIND_ARGS=()
if [[ -n "${BIND}" ]]; then
  SING_BIND_ARGS=(--bind "${BIND}")
fi

# Stage output on node-local disk when available, then copy to OUTDIR (NFS-friendly).
STAGE="${TMPDIR:-/tmp}/kmat_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${STAGE}"
STAGE_KSET="${STAGE}/${STEM}.kset"

singularity exec "${SING_BIND_ARGS[@]}" "${KMAT_SIF}" \
  kmat --profile "${KMAT_PROFILE}" --threads "${SLURM_CPUS_PER_TASK}" \
  count -i "${FASTQ}" -s "${KMER_SIZE}" --ci "${CI}" -o "${STAGE_KSET}"

cp -f "${STAGE_KSET}" "${KSET}"
rm -rf "${STAGE}"

echo "[$(date -Is)] done ${KSET}"
