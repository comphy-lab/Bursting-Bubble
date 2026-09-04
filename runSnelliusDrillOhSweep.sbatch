#!/bin/bash
#SBATCH --job-name=sbbDrill
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=72:00:00
#SBATCH --partition=genoa
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vatsal.sanjay@comphy-lab.org
#SBATCH --output=slurm-drill-c%x-%j.out
#SBATCH --error=slurm-drill-c%x-%j.err

# =============================================================================
# Bursting Bubble - drill (feature-tracking AMR) single-case runner on Snellius.
#
# One Slurm job = one Oh-sweep case, run FROM SCRATCH with the full production
# regularisation stack (focus cap, gas-wisp removal, arm/fire inception latch,
# regional post-latch ceiling, terminal relaxation) — see DRILL-RESOLUTION.md
# in the repo root for the mechanism and the validated knob values this script
# hard-codes below.
#
# The drill solver is OpenMP-only (distance.h init is MPI-incompatible), so
# this is a single-node, single-task, multi-threaded job — same pattern as
# the ElasticPinchOff / VE-drop-impact OpenMP jobs elsewhere on this account
# (see e.g. /projects/0/nctt0620/vatsal/2026-06-27-Elastic-Pinch-Off/runSnellius_c*.sbatch).
#
# Usage (one job per Oh value in the sweep):
#   sbatch --job-name=sbb2001 \
#          --export=ALL,CASE_NO=2001,OH=0.0216,MAXLEVEL=13 \
#          runSnelliusDrillOhSweep.sbatch
#
# Required env (via --export): CASE_NO, OH
# Optional env (defaults shown): MAXLEVEL=13, BOND=1e-3, TMAX=1.5,
#   DRILL_MAXLEVEL_START=8, DRILL_MAXLEVEL_FOCUS=<MAXLEVEL-1>,
#   DRILL_RELAX_LEVEL=10, DRILL_REMOVE_GAS_SIZE=3, KE_STOP_MAX=1e4,
#   TIP_METRICS_LOG=0, DRILL_HOLD_MAX_UNTIL_TIP_PINCH=0,
#   DRILL_STOP_AT_TIP_PINCH=0, DRILL_TSNAP_STAGES=1
# =============================================================================

set -euo pipefail

: "${CASE_NO:?CASE_NO must be set (--export=ALL,CASE_NO=...)}"
: "${OH:?OH must be set (--export=ALL,OH=...)}"
MAXLEVEL="${MAXLEVEL:-13}"
BOND="${BOND:-1e-3}"
TMAX="${TMAX:-1.5}"
DRILL_MAXLEVEL_START="${DRILL_MAXLEVEL_START:-8}"
DRILL_MAXLEVEL_FOCUS="${DRILL_MAXLEVEL_FOCUS:-$((MAXLEVEL - 1))}"
DRILL_RELAX_LEVEL="${DRILL_RELAX_LEVEL:-10}"
DRILL_REMOVE_GAS_SIZE="${DRILL_REMOVE_GAS_SIZE:-3}"
DRILL_TSNAP_STAGES="${DRILL_TSNAP_STAGES:-1}"
TIP_METRICS_LOG="${TIP_METRICS_LOG:-0}"
DRILL_HOLD_MAX_UNTIL_TIP_PINCH="${DRILL_HOLD_MAX_UNTIL_TIP_PINCH:-0}"
DRILL_STOP_AT_TIP_PINCH="${DRILL_STOP_AT_TIP_PINCH:-0}"
KE_STOP_MAX="${KE_STOP_MAX:-1e4}"
[[ "${CASE_NO}" =~ ^[0-9]{4}$ ]] || { echo "CASE_NO must be four digits" >&2; exit 2; }
[[ "${MAXLEVEL}" =~ ^[0-9]+$ ]] || { echo "MAXLEVEL must be a nonnegative integer" >&2; exit 2; }
[[ "${DRILL_MAXLEVEL_START}" =~ ^[0-9]+$ ]] || { echo "DRILL_MAXLEVEL_START must be a nonnegative integer" >&2; exit 2; }
[[ "${DRILL_MAXLEVEL_FOCUS}" =~ ^-?[0-9]+$ ]] || { echo "DRILL_MAXLEVEL_FOCUS must be an integer" >&2; exit 2; }
for flag in DRILL_TSNAP_STAGES TIP_METRICS_LOG DRILL_HOLD_MAX_UNTIL_TIP_PINCH DRILL_STOP_AT_TIP_PINCH; do
  [[ "${!flag}" == 0 || "${!flag}" == 1 ]] || { echo "${flag} must be 0 or 1" >&2; exit 2; }
done
# staged-tsnap floor: 2^(start-MAXlevel), so the fast jet window near
# inception is never under-sampled at the target resolution
DRILL_TSNAP_MIN_FACTOR=$(python3 -c "print(2**(${DRILL_MAXLEVEL_START}-${MAXLEVEL}))")

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
THREADS="${SLURM_CPUS_PER_TASK:-24}"
cd "$SCRIPT_DIR"

echo "============================================="
echo "Bursting Bubble - drill - Snellius Oh-sweep case ${CASE_NO}"
echo "============================================="
echo "Started:   $(date)"
echo "Node:      $(hostname)"
echo "Job ID:    ${SLURM_JOB_ID:-local}"
echo "Threads:   ${THREADS}"
echo "Oh=${OH} Bond=${BOND} MAXlevel=${MAXLEVEL} focus-cap=${DRILL_MAXLEVEL_FOCUS} tmax=${TMAX}"
echo ""

echo "Loading modules..."
module purge
module load 2024
module load OpenMPI/5.0.3-GCC-13.3.0
echo ""

# Pin OpenMP threads to cores: shared genoa nodes, keep threads from
# migrating across NUMA domains during a multi-hour run.
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# Reuse the account's shared Basilisk install (already ref-locked and
# verified on this cluster) rather than installing a fresh per-project copy.
export BASILISK="${HOME}/basilisk/src"
if [[ ! -x "${BASILISK}/qcc" ]]; then
  echo "ERROR: qcc not found at ${BASILISK}/qcc — Basilisk not installed at \$HOME/basilisk." >&2
  exit 1
fi
export PATH="${BASILISK}:${PATH}"

# --- Case directory + params ---------------------------------------------
CASE_DIR="${SCRIPT_DIR}/simulationCases/${CASE_NO}"
if [[ -e "${CASE_DIR}" ]]; then
  [[ -d "${CASE_DIR}" ]] || {
    echo "ERROR: case path exists and is not a directory: ${CASE_DIR}" >&2
    exit 2
  }
  [[ -z "$(find "${CASE_DIR}" -mindepth 1 -maxdepth 1 -print -quit)" ]] || {
    echo "ERROR: case directory is not empty; refusing an implicit restart: ${CASE_DIR}" >&2
    exit 2
  }
else
  mkdir "${CASE_DIR}"
fi
cp "${SCRIPT_DIR}/simulationCases/burstingBubble-drillResolution.c" "${CASE_DIR}/"

cat > "${CASE_DIR}/case.params" <<PARAMS
CaseNo=${CASE_NO}
Oh=${OH}
Bond=${BOND}
OhRatio=2e-2
MAXlevel=${MAXLEVEL}
MINlevel=4
init_grid_level=5
fErr=1e-3
VelErr=1e-3
KErr=1e-6
CFL=0.1
dtmax=1e-2
TOLERANCE=1e-4
tsnap=1e-2
tmax=${TMAX}
zWall=4.0
drillAMR=1
drillMaxlevelStart=${DRILL_MAXLEVEL_START}
drillNcellsK=5
drillNcellsJet=5
drillRelaxLevel=${DRILL_RELAX_LEVEL}
drillTsnapStages=${DRILL_TSNAP_STAGES}
drillTsnapMinFactor=${DRILL_TSNAP_MIN_FACTOR}
drillMaxlevelFocus=${DRILL_MAXLEVEL_FOCUS}
drillRemoveGasSize=${DRILL_REMOVE_GAS_SIZE}
tipMetricsLog=${TIP_METRICS_LOG}
drillHoldMaxUntilTipPinch=${DRILL_HOLD_MAX_UNTIL_TIP_PINCH}
drillStopAtTipPinch=${DRILL_STOP_AT_TIP_PINCH}
keStopMax=${KE_STOP_MAX}
PARAMS
echo "Wrote ${CASE_DIR}/case.params:"
cat "${CASE_DIR}/case.params"
echo ""

# --- Compile ----------------------------------------------------------------
cd "${CASE_DIR}"
echo "Compiling burstingBubble-drillResolution.c ..."
qcc -O2 -Wall -disable-dimensions -fopenmp -I../../src-local \
    burstingBubble-drillResolution.c -o drill -lm

# --- Run (from scratch: distance.h init, so this must stay single-node) ----
echo "Launching case ${CASE_NO} (Oh=${OH}, OpenMP x${THREADS}) ..."
rc=0
OMP_NUM_THREADS="${THREADS}" ./drill case.params || rc=$?

echo ""
echo "============================================="
echo "Finished:  $(date)  rc=${rc}"
echo "Output:    ${CASE_DIR}/"
echo "Log:       ${CASE_DIR}/log"
echo "============================================="
exit "$rc"
