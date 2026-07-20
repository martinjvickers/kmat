#!/usr/bin/env bash
# Submit count array + dependent build for a FASTQ accession list.
#
# Usage:
#   export KMAT_SIF=/path/to/kmat.sif
#   export ACCESSION_LIST=/path/to/my_100_fastqs.txt
#   export OUTDIR=/scratch/$USER/kmat_run1
#   export BIND=/project:/project,/scratch:/scratch   # adjust to your site
#   export CI=2
#   export KMER_SIZE=31
#   ./hpc/submit_panel.sh
#
# ACCESSION_LIST format: one FASTQ/FASTA path per line (see accession_list.example.txt).

set -euo pipefail

: "${KMAT_SIF:?}"
: "${ACCESSION_LIST:?}"
: "${OUTDIR:?}"

BIND="${BIND:-}"
CI="${CI:-2}"
KMER_SIZE="${KMER_SIZE:-31}"
KMAT_PROFILE="${KMAT_PROFILE:-hpc}"
ARRAY_THROTTLE="${ARRAY_THROTTLE:-50}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "${OUTDIR}/ksets" logs

N=$(grep -vE '^\s*(#|$)' "${ACCESSION_LIST}" | wc -l | tr -d ' ')
if [[ "${N}" -lt 1 ]]; then
  echo "ERROR: no accessions in ${ACCESSION_LIST}" >&2
  exit 1
fi

echo "Submitting count array for ${N} accessions → ${OUTDIR}/ksets/"
EXPORT="ALL,KMAT_SIF=${KMAT_SIF},ACCESSION_LIST=${ACCESSION_LIST},OUTDIR=${OUTDIR},CI=${CI},KMER_SIZE=${KMER_SIZE},KMAT_PROFILE=${KMAT_PROFILE},BIND=${BIND}"

COUNT_JOB=$(sbatch --parsable \
  --array="1-${N}%${ARRAY_THROTTLE}" \
  --export="${EXPORT}" \
  "${SCRIPT_DIR}/slurm_count_array.sh")
echo "count job: ${COUNT_JOB}"

BUILD_JOB=$(sbatch --parsable \
  --dependency="afterok:${COUNT_JOB}" \
  --export="${EXPORT}" \
  "${SCRIPT_DIR}/slurm_build.sh")
echo "build job: ${BUILD_JOB} (afterok:${COUNT_JOB})"

echo
echo "Monitor:"
echo "  squeue -u \$USER"
echo "  tail -f logs/kmat_count_${COUNT_JOB}_1.out"
echo "  ls ${OUTDIR}/ksets | wc -l"
echo "Matrix will be: ${OUTDIR}/panel.kmat"
