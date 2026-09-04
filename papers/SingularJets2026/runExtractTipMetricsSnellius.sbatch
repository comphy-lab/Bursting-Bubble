#!/bin/bash
#SBATCH --job-name=sbbtip
#SBATCH --account=nctth181
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=04:00:00
#SBATCH --partition=genoa
#SBATCH --output=slurm-tip-metrics-%x-%j.out
#SBATCH --error=slurm-tip-metrics-%x-%j.err

# Extract tip curvature and velocity from one case. Production mode uses every
# core on one Genoa node; bounded smoke mode may share a node but requires an
# explicit MAX_FRAMES cap.
set -euo pipefail

: "${CASEDIR:?CASEDIR must be set}"
: "${OUT:?OUT must be set}"
: "${SOURCE_ROOT:?SOURCE_ROOT must name the immutable source checkout}"
: "${OH:?OH must be set}"
: "${BOND:?BOND must be set}"
: "${PRE_LEVEL:?PRE_LEVEL must be set}"
: "${POST_LEVEL:?POST_LEVEL must be set}"
: "${T0:?T0 must be set}"
: "${T0_PROTOCOL:?T0_PROTOCOL must be set}"
: "${LINEAGE:?LINEAGE must be set}"

RUN_ROOT="${RUN_ROOT:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
RUN_REAL="$(realpath "${RUN_ROOT}")"
SOURCE_REAL="$(realpath "${SOURCE_ROOT}")"
CASE_REAL="$(realpath "${CASEDIR}")"
OUT_REAL="$(realpath -m "${OUT}")"
[[ "${SOURCE_REAL}" == "${RUN_REAL}"/postprocess-source/* ]] || {
    echo "SOURCE_ROOT must be below ${RUN_REAL}/postprocess-source" >&2
    exit 2
}
[[ "${CASE_REAL}" == "${RUN_REAL}"/simulationCases/* ]] || {
    echo "CASEDIR must be below ${RUN_REAL}/simulationCases" >&2
    exit 2
}
[[ "${OUT_REAL}" == "${RUN_REAL}"/postprocessed/* ]] || {
    echo "OUT must be below ${RUN_REAL}/postprocessed" >&2
    exit 2
}

THREADS="${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is unset}"
MODE="${MODE:-production}"
FORCE="${FORCE:-0}"
BRIDGED="${BRIDGED:-0}"
[[ "${FORCE}" == 0 || "${FORCE}" == 1 ]] || { echo "FORCE must be 0 or 1" >&2; exit 2; }
[[ "${BRIDGED}" == 0 || "${BRIDGED}" == 1 ]] || { echo "BRIDGED must be 0 or 1" >&2; exit 2; }
case "${MODE}" in
    production)
        [[ "${THREADS}" -eq 192 ]] || {
            echo "production extraction requires 192 CPUs" >&2
            exit 2
        }
        ;;
    smoke)
        [[ "${THREADS}" -gt 0 && "${THREADS}" -le 24 ]] || {
            echo "smoke extraction requires 1-24 CPUs" >&2
            exit 2
        }
        [[ "${MAX_FRAMES:-}" =~ ^[1-9][0-9]*$ ]] || {
            echo "smoke extraction requires positive MAX_FRAMES" >&2
            exit 2
        }
        ;;
    *)
        echo "MODE must be production or smoke" >&2
        exit 2
        ;;
esac

module purge
module load 2024
module load Python/3.12.3-GCCcore-13.3.0

export BASILISK="${BASILISK:-${HOME}/basilisk/src}"
export PATH="${BASILISK}:${PATH}"
test -x "${BASILISK}/qcc"

EXTRACTOR="${SOURCE_REAL}/papers/SingularJets2026/tip-curvature/extract_tip_metrics.py"
TIP_SOURCE="${SOURCE_REAL}/postProcess/getTipMetrics.c"
test -f "${EXTRACTOR}"
test -f "${TIP_SOURCE}"

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/sbb-tip-${SLURM_JOB_ID}.XXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT
cp "${TIP_SOURCE}" "${WORK_DIR}/getTipMetrics.c"
pushd "${WORK_DIR}" >/dev/null
qcc -O2 -Wall -disable-dimensions getTipMetrics.c -o getTipMetrics -lm
popd >/dev/null

ARGS=(
    --case "${CASE_REAL}"
    --out "${OUT_REAL}"
    --gettip "${WORK_DIR}/getTipMetrics"
    --cpus "${THREADS}"
    --oh "${OH}"
    --bond "${BOND}"
    --pre-level "${PRE_LEVEL}"
    --post-level "${POST_LEVEL}"
    --t0 "${T0}"
    --t0-protocol "${T0_PROTOCOL}"
    --lineage "${LINEAGE}"
    --skip-video
)
[[ -n "${TIME_MIN:-}" ]] && ARGS+=(--time-min "${TIME_MIN}")
[[ -n "${TIME_MAX:-}" ]] && ARGS+=(--time-max "${TIME_MAX}")
[[ -n "${MAX_FRAMES:-}" ]] && ARGS+=(--max-frames "${MAX_FRAMES}")
[[ "${FORCE}" == 1 ]] && ARGS+=(--force)
[[ "${BRIDGED}" == 1 ]] && ARGS+=(--bridged)

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
srun --exclusive --cpu-bind=cores --ntasks=1 --cpus-per-task="${THREADS}" \
    python3 "${EXTRACTOR}" "${ARGS[@]}"

test -s "${OUT_REAL}/tip_metrics.csv"
test -s "${OUT_REAL}/manifest.json"
echo "TIP_METRICS_COMPLETE out=${OUT_REAL}"
