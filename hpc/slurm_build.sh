#!/bin/bash
#SBATCH --job-name=kmat-build
#SBATCH --output=logs/kmat_build_%j.out
#SBATCH --error=logs/kmat_build_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#
# Build a v2 PA matrix from .kset files produced by slurm_count_array.sh.
#
# Required:
#   KMAT_SIF          path to kmat.sif
#   ACCESSION_LIST    same FASTQ list used for counting (order = matrix columns)
#   OUTDIR            directory containing ksets/ from the count array
# Optional:
#   KMER_SIZE         default 31
#   BIND              singularity bind specs
#   KMAT_PROFILE      default hpc
#   MATRIX_NAME       default panel.kmat
#
# Submit after the count array finishes, e.g.:
#   sbatch --dependency=afterok:<count_jobid> \
#     --export=ALL,KMAT_SIF=...,ACCESSION_LIST=...,OUTDIR=...,BIND=/data:/data \
#     hpc/slurm_build.sh

set -euo pipefail

: "${KMAT_SIF:?set KMAT_SIF}"
: "${ACCESSION_LIST:?set ACCESSION_LIST}"
: "${OUTDIR:?set OUTDIR}"

KMER_SIZE="${KMER_SIZE:-31}"
KMAT_PROFILE="${KMAT_PROFILE:-hpc}"
BIND="${BIND:-}"
MATRIX_NAME="${MATRIX_NAME:-panel.kmat}"

mkdir -p "${OUTDIR}" logs

KSET_LIST="${OUTDIR}/kset_list.txt"
: > "${KSET_LIST}"

mapfile -t LINES < <(grep -vE '^\s*(#|$)' "${ACCESSION_LIST}")
if [[ "${#LINES[@]}" -eq 0 ]]; then
  echo "ERROR: empty accession list" >&2
  exit 1
fi

missing=0
for FASTQ in "${LINES[@]}"; do
  BASE="$(basename "${FASTQ}")"
  STEM="${BASE}"
  for ext in .fastq.gz .fq.gz .fasta.gz .fa.gz .fna.gz .fastq .fq .fasta .fa .fna .gz; do
    case "${STEM}" in
      *"${ext}") STEM="${STEM%${ext}}"; break ;;
    esac
  done
  KSET="${OUTDIR}/ksets/${STEM}.kset"
  if [[ ! -f "${KSET}" ]]; then
    echo "MISSING: ${KSET} (from ${FASTQ})" >&2
    missing=$((missing + 1))
    continue
  fi
  echo "${KSET}" >> "${KSET_LIST}"
done

if [[ "${missing}" -ne 0 ]]; then
  echo "ERROR: ${missing} .kset file(s) missing; not building" >&2
  exit 1
fi

MATRIX="${OUTDIR}/${MATRIX_NAME}"
echo "[$(date -Is)] building ${MATRIX} from ${#LINES[@]} presence sets"

SING_BIND_ARGS=()
if [[ -n "${BIND}" ]]; then
  SING_BIND_ARGS=(--bind "${BIND}")
fi

STAGE="${TMPDIR:-/tmp}/kmat_build_${SLURM_JOB_ID}"
mkdir -p "${STAGE}"
STAGE_MATRIX="${STAGE}/${MATRIX_NAME}"
STAGE_LIST="${STAGE}/kset_list.txt"
cp "${KSET_LIST}" "${STAGE_LIST}"

# Bind OUTDIR so singularity can read .kset paths (absolute host paths).
# Prefer staging the list with paths that remain valid inside the container via BIND.
singularity exec "${SING_BIND_ARGS[@]}" "${KMAT_SIF}" \
  kmat --profile "${KMAT_PROFILE}" --threads "${SLURM_CPUS_PER_TASK}" \
  build -k "${KSET_LIST}" -s "${KMER_SIZE}" -o "${STAGE_MATRIX}"

cp -f "${STAGE_MATRIX}" "${MATRIX}"
rm -rf "${STAGE}"

echo "[$(date -Is)] wrote ${MATRIX}"
echo "Next (optional):"
echo "  singularity exec \$BIND_ARGS ${KMAT_SIF} kmat validate -i ${MATRIX} -k ${KSET_LIST}"
echo "  singularity exec \$BIND_ARGS ${KMAT_SIF} kmat pop -i ${MATRIX} -k ${KSET_LIST} -o ${OUTDIR}/pop.tsv"
