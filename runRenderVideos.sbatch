#!/bin/bash
#SBATCH --job-name=sbbvid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=04:00:00
#SBATCH --partition=genoa
#SBATCH --output=slurm-vid-%x-%j.out
# =============================================================================
# Render the two post-process videos for ONE case ON the compute node, from the
# case's local snapshot dumps — no data transfer off the cluster (only the small
# MP4s leave). Companion to runSnelliusDrillOhSweep.sbatch.
#
#   Video 1 (vel+dissipation): postProcess/Video.py   -> velocity + strain rate
#   Video 2 (interface+mesh):  postProcess/render_pair_video.py
#                              -> interface + jet-base marker | interface + AMR
#
# Both scripts shell out to the compiled helpers in postProcess/:
#   getData, getView2D           (compile clean: qcc -O2 -disable-dimensions ...)
#   getFacet -> use getFacetFull (bare getFacet.c segfaults on two-phase axi
#               dumps; build getFacetFull.c and install it as ./getFacet)
#   getBase                      (bare getBase.c is fine on same-version dumps)
# getView2D uses basilisk's software framebuffer (libfb_tiny) — no GPU needed.
#
# Usage (one job per case):
#   sbatch --job-name=sbbvid5001 \
#          --export=ALL,CASE=5001,CASEDIR=/abs/path/simulationCases/5001,TAG=L15 \
#          runRenderVideos.sbatch
#
# Required env (via --export): CASE, CASEDIR, TAG (a label -> <CASE>_<TAG>_*.mp4)
# Optional env: POSTPROCESS (default <submit-dir>/postProcess), VOUT (default
#   <submit-dir>/videos), FRAMING "ZMIN ZMAX RMAX FPS" (default cavity+jet).
# =============================================================================
set -uo pipefail
: "${CASE:?CASE must be set}"; : "${CASEDIR:?CASEDIR must be set}"; : "${TAG:?TAG must be set}"

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
PP="${POSTPROCESS:-${SUBMIT_DIR}/postProcess}"
VOUT="${VOUT:-${SUBMIT_DIR}/videos}"
mkdir -p "${VOUT}"

# --- python (numpy+matplotlib+Pillow) + ffmpeg + Basilisk. Module names below
#     are for Snellius/RHEL9 2024; adjust for your site. ------------------------
module load 2024 matplotlib/3.9.2-gfbf-2024a Pillow/10.4.0-GCCcore-13.3.0 \
            FFmpeg/7.0.2-GCCcore-13.3.0 2>/dev/null || true
export BASILISK="${BASILISK:-$HOME/basilisk/src}"; export PATH="${BASILISK}:${PATH}"

THREADS="${SLURM_CPUS_PER_TASK:-24}"
read -r ZMIN ZMAX RMAX FPS <<< "${FRAMING:--2.05 1.5 1.1 18}"
cd "${PP}"

# Video 1: velocity + strain-rate (dissipation). Video.py writes <CASEDIR>/<CASE>.mp4
python3 Video.py --caseToProcess "${CASEDIR}" --CPUs "${THREADS}" \
        --ZMIN "${ZMIN}" --ZMAX "${ZMAX}" --RMAX "${RMAX}" \
        --folderToSave "${VOUT}/frames_${CASE}_vel"
[ -f "${CASEDIR}/${CASE}.mp4" ] && cp "${CASEDIR}/${CASE}.mp4" \
        "${VOUT}/${CASE}_${TAG}_vel_dissipation.mp4"

# Video 2: interface + jet-base marker | interface + AMR cells (all snapshots)
python3 render_pair_video.py "${CASEDIR}" 10 \
        "${VOUT}/${CASE}_${TAG}_interface_marker_AMR.mp4" "${ZMIN}" "${ZMAX}" "${RMAX}" "${FPS}"

echo "DONE ${CASE}" \
     "v1=$([ -f ${VOUT}/${CASE}_${TAG}_vel_dissipation.mp4 ] && echo ok || echo MISS)" \
     "v2=$([ -f ${VOUT}/${CASE}_${TAG}_interface_marker_AMR.mp4 ] && echo ok || echo MISS)"
