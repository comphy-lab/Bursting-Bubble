#!/bin/bash
#SBATCH --job-name=sbbfacets
#SBATCH --account=nctth181
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=04:00:00
#SBATCH --partition=genoa
#SBATCH --exclusive
#SBATCH --output=slurm-facets-%x-%j.out
#SBATCH --error=slurm-facets-%x-%j.err
# =============================================================================
# Full-node interface-facet extraction for one completed simulation case.
#
# Required environment: CASEDIR, OUT, T0.
# Optional: PRE_LO, RMAX (0.6), MAX_LOG_GAP (1e-3), MAX_FRAMES, FORCE (0|1).
# Submit from the Bursting-Bubble run root so CASEDIR remains under its
# simulationCases/ tree and OUT remains under postprocessed/.
# =============================================================================
set -euo pipefail

: "${CASEDIR:?CASEDIR must be set}"
: "${OUT:?OUT must be set}"
: "${T0:?T0 must be set}"

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
SUBMIT_REAL="$(realpath "${SUBMIT_DIR}")"
CASE_REAL="$(realpath "${CASEDIR}")"
OUT_REAL="$(realpath -m "${OUT}")"
[[ "${CASE_REAL}" == "${SUBMIT_REAL}"/simulationCases/* ]] || {
    echo "CASEDIR must be below ${SUBMIT_REAL}/simulationCases" >&2
    exit 2
}
[[ "${OUT_REAL}" == "${SUBMIT_REAL}"/postprocessed/* ]] || {
    echo "OUT must be below ${SUBMIT_REAL}/postprocessed" >&2
    exit 2
}

THREADS="${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is unset}"
[[ "${THREADS}" -eq 192 ]] || { echo "Facet extraction requires 192 CPUs" >&2; exit 2; }
RMAX="${RMAX:-0.6}"
MAX_LOG_GAP="${MAX_LOG_GAP:-1e-3}"
FORCE="${FORCE:-0}"
[[ "${FORCE}" == 0 || "${FORCE}" == 1 ]] || { echo "FORCE must be 0 or 1" >&2; exit 2; }

module purge
module load 2024
module load Python/3.12.3-GCCcore-13.3.0

export BASILISK="${BASILISK:-${HOME}/basilisk/src}"
export PATH="${BASILISK}:${PATH}"
test -x "${BASILISK}/qcc"

EXTRACTOR="${SUBMIT_REAL}/papers/SingularJets2026/facet-collapse-figure3/extract_full.py"
FACET_SOURCE="${SUBMIT_REAL}/postProcess/getFacetFull.c"
test -f "${EXTRACTOR}"
test -f "${FACET_SOURCE}"

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/sbb-facets-${SLURM_JOB_ID}.XXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT
cp "${FACET_SOURCE}" "${WORK_DIR}/getFacetFull.c"
pushd "${WORK_DIR}" >/dev/null
qcc -O2 -Wall -disable-dimensions getFacetFull.c -o getFacet -lm
popd >/dev/null

ARGS=(
    --case "${CASE_REAL}"
    --out "${OUT_REAL}"
    --getfacet "${WORK_DIR}/getFacet"
    --t0 "${T0}"
    --rmax "${RMAX}"
    --max-log-gap "${MAX_LOG_GAP}"
    --cpus "${THREADS}"
)
[[ -n "${PRE_LO:-}" ]] && ARGS+=(--pre-lo "${PRE_LO}")
[[ -n "${MAX_FRAMES:-}" ]] && ARGS+=(--max-frames "${MAX_FRAMES}")
[[ "${FORCE}" == 1 ]] && ARGS+=(--force)

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
srun --exclusive --cpu-bind=cores --ntasks=1 --cpus-per-task="${THREADS}" \
    python3 "${EXTRACTOR}" "${ARGS[@]}"

test -s "${OUT_REAL}/index_pre.txt"
test -s "${OUT_REAL}/index.txt"
test -s "${OUT_REAL}/manifest.json"
echo "FACETS_COMPLETE out=${OUT_REAL}"
