#!/bin/bash
#SBATCH --job-name=sbb-tip-sweep
#SBATCH --account=nctth181
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=5-00:00:00
#SBATCH --partition=genoa
#SBATCH --array=0-7%8
#SBATCH --output=slurm-tip-sweep-%A_%a.out
#SBATCH --error=slurm-tip-sweep-%A_%a.err

# One logarithmic-Oh short-window curvature case per array task. Submit tasks
# 0 and 6 first as endpoint smokes. Their post-run acceptance gate exits
# nonzero unless the requested tip level is reached and a resolved, interior,
# rebounding pre-pinch minimum is present; only then release the remaining
# array indices through an afterok dependency.
set -euo pipefail

ROOT="$(realpath "${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is unset}")"
SOURCE_ROOT="${SOURCE_ROOT:-${ROOT}}"
SOURCE_REAL="$(realpath "${SOURCE_ROOT}")"
[[ "${SOURCE_REAL}" == "${ROOT}" || "${SOURCE_REAL}" == "${ROOT}"/postprocess-source/* ]] || {
    echo "SOURCE_ROOT must be the run root or an immutable checkout below postprocess-source" >&2
    exit 2
}
TABLE="${SWEEP_TABLE:-${SOURCE_REAL}/papers/SingularJets2026/tip-curvature/logspace-sweep.csv}"
RUNNER="${SOURCE_REAL}/runSnelliusDrillOhSweep.sbatch"
VALIDATOR="${SOURCE_REAL}/papers/SingularJets2026/tip-curvature/validate_online_tip_metrics.py"
test -f "${TABLE}"
test -f "${RUNNER}"
test -f "${VALIDATOR}"

ARRAY_INDEX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is unset}"
[[ "${ARRAY_INDEX}" =~ ^[0-7]$ ]] || { echo "array index must be 0..7" >&2; exit 2; }
ROW="$(awk -F, -v wanted="${ARRAY_INDEX}" '
  NR > 1 && $1 == wanted {print; found++}
  END {if (found != 1) exit 3}
' "${TABLE}")" || { echo "sweep table must contain one row for index ${ARRAY_INDEX}" >&2; exit 2; }
IFS=, read -r TABLE_INDEX CASE_NO OH POST_LEVEL PRE_LEVEL <<< "${ROW}"
[[ "${TABLE_INDEX}" == "${ARRAY_INDEX}" ]] || exit 2
[[ "${CASE_NO}" =~ ^[0-9]{4}$ ]] || { echo "invalid case number in sweep table" >&2; exit 2; }
[[ "${OH}" =~ ^0[.][0-9]+$ ]] || { echo "invalid Oh in sweep table" >&2; exit 2; }
[[ "${POST_LEVEL}" =~ ^(18|19)$ ]] || { echo "invalid post level in sweep table" >&2; exit 2; }
[[ "${PRE_LEVEL}" == 12 ]] || { echo "pre level must be 12" >&2; exit 2; }

export CASE_NO OH MAXLEVEL="${POST_LEVEL}" BOND=0 TMAX=0.60
export DRILL_MAXLEVEL_START=8 DRILL_MAXLEVEL_FOCUS="${PRE_LEVEL}"
export DRILL_RELAX_LEVEL=10 DRILL_REMOVE_GAS_SIZE=3 KE_STOP_MAX=1e6
export DRILL_TSNAP_STAGES=0 TIP_METRICS_LOG=1
export DRILL_HOLD_MAX_UNTIL_TIP_PINCH=1 DRILL_STOP_AT_TIP_PINCH=1

cd "${ROOT}"
bash "${RUNNER}"

module load Python/3.12.3-GCCcore-13.3.0
python3 "${VALIDATOR}" \
    --log "${ROOT}/simulationCases/${CASE_NO}/tip_metrics.log" \
    --target-level "${POST_LEVEL}" \
    --min-valid-rows 50 \
    --min-cells 4 \
    --rebound-factor 1.1 \
    --rebound-rows 20 \
    --output-json "${ROOT}/simulationCases/${CASE_NO}/tip_metrics_gate.json"
