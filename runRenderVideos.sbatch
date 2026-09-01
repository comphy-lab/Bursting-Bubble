#!/bin/bash
#SBATCH --job-name=sbbvidframes
#SBATCH --account=nctth181
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --partition=genoa
#SBATCH --exclusive
#SBATCH --output=slurm-video-frames-%x-%j.out
#SBATCH --error=slurm-video-frames-%x-%j.err
# =============================================================================
# Render both expensive frame sets for one completed case on one full Snellius
# genoa node. Encoding is deliberately split into runEncodeVideos.sbatch so the
# 192-core allocation never spends its useful life running a serial ffmpeg.
#
# Required environment: CASE, CASEDIR, TAG, VOUT.
# Optional: FRAMING="ZMIN ZMAX RMAX FPS", GRIDS_PER_R (default 256).
# =============================================================================
set -euo pipefail

: "${CASE:?CASE must be set}"
: "${CASEDIR:?CASEDIR must be set}"
: "${TAG:?TAG must be set}"
: "${VOUT:?VOUT must be set}"

[[ "${CASE}" =~ ^[0-9]{4}$ ]] || { echo "CASE must be four digits" >&2; exit 2; }
[[ "${TAG}" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "TAG contains unsafe characters" >&2; exit 2; }

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
SUBMIT_REAL="$(realpath "${SUBMIT_DIR}")"
CASE_REAL="$(realpath "${CASEDIR}")"
VOUT_REAL="$(realpath -m "${VOUT}")"
[[ "${CASE_REAL}" == "${SUBMIT_REAL}"/simulationCases/"${CASE}" ]] || {
    echo "CASEDIR must be simulationCases/${CASE} below the submit root" >&2
    exit 2
}
[[ "${VOUT_REAL}" == "${SUBMIT_REAL}"/postprocessed/* ]] || {
    echo "VOUT must be below ${SUBMIT_REAL}/postprocessed" >&2
    exit 2
}

THREADS="${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK is unset}"
[[ "${THREADS}" -eq 192 ]] || { echo "Frame job requires all 192 node CPUs" >&2; exit 2; }
GRIDS_PER_R="${GRIDS_PER_R:-256}"
[[ "${GRIDS_PER_R}" =~ ^[0-9]+$ ]] && [[ "${GRIDS_PER_R}" -gt 0 ]] || {
    echo "GRIDS_PER_R must be positive" >&2
    exit 2
}
read -r ZMIN ZMAX RMAX FPS <<< "${FRAMING:--2.05 1.5 1.1 18}"

module purge
module load 2024
module load Python/3.12.3-GCCcore-13.3.0
module load SciPy-bundle/2024.05-gfbf-2024a
module load matplotlib/3.9.2-gfbf-2024a
module load Pillow/10.4.0-GCCcore-13.3.0

export BASILISK="${BASILISK:-${HOME}/basilisk/src}"
export PATH="${BASILISK}:${PATH}"
test -x "${BASILISK}/qcc"
python3 -c 'import matplotlib, numpy, PIL'

PP_SOURCE="${SUBMIT_REAL}/postProcess"
for source in Video.py render_pair_video.py getFacetFull.c getData.c getBase.c getView2D.c; do
    test -f "${PP_SOURCE}/${source}"
done

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/sbb-video-${SLURM_JOB_ID}.XXXXXX")"
trap 'rm -rf -- "${WORK_DIR}"' EXIT
cp "${PP_SOURCE}/Video.py" "${PP_SOURCE}/render_pair_video.py" "${WORK_DIR}/"

qcc -O2 -Wall -disable-dimensions "${PP_SOURCE}/getFacetFull.c" -o "${WORK_DIR}/getFacet" -lm
qcc -O2 -Wall -disable-dimensions "${PP_SOURCE}/getData.c" -o "${WORK_DIR}/getData" -lm
qcc -O2 -Wall -disable-dimensions "${PP_SOURCE}/getBase.c" -o "${WORK_DIR}/getBase" -lm
qcc -O2 -Wall -disable-dimensions "${PP_SOURCE}/getView2D.c" -o "${WORK_DIR}/getView2D" \
    -L"${BASILISK}/gl" -lglutils -lfb_tiny -lm

SNAPSHOT_COUNT="$(find "${CASE_REAL}/intermediate" -maxdepth 1 -type f -name 'snapshot-*' -printf . | wc -c)"
ZERO_COUNT="$(find "${CASE_REAL}/intermediate" -maxdepth 1 -type f -name 'snapshot-*' -size 0 -printf . | wc -c)"
[[ "${SNAPSHOT_COUNT}" -gt 0 ]] || { echo "No snapshots found" >&2; exit 3; }
[[ "${ZERO_COUNT}" -eq 0 ]] || { echo "Zero-byte snapshots found" >&2; exit 3; }
test -s "${CASE_REAL}/log"
test -s "${CASE_REAL}/case.params"

mkdir -p "${VOUT_REAL}"
VEL_FRAMES="${VOUT_REAL}/frames_${CASE}_vel"
PAIR_FRAMES="${VOUT_REAL}/frames_${CASE}_pair"
CONFIG_FILE="${VOUT_REAL}/${CASE}_${TAG}_render-config.txt"
CONFIG_VALUE="case=${CASE} tag=${TAG} snapshots=${SNAPSHOT_COUNT} framing=${ZMIN},${ZMAX},${RMAX},${FPS} grids_per_r=${GRIDS_PER_R}"
if [[ -f "${CONFIG_FILE}" ]]; then
    [[ "$(<"${CONFIG_FILE}")" == "${CONFIG_VALUE}" ]] || {
        echo "Existing render configuration differs; use a new VOUT" >&2
        exit 3
    }
else
    printf '%s\n' "${CONFIG_VALUE}" > "${CONFIG_FILE}"
fi

# Four-frame serial/parallel equivalence catches helper incompatibility,
# shared-scratch races and non-deterministic output before the full-node pass.
PILOT_ONE="${WORK_DIR}/pilot-one"
PILOT_FOUR="${WORK_DIR}/pilot-four"
mkdir -p "${PILOT_ONE}/vel" "${PILOT_ONE}/pair" "${PILOT_FOUR}/vel" "${PILOT_FOUR}/pair"
python3 "${WORK_DIR}/Video.py" --caseToProcess "${CASE_REAL}" --CPUs 1 --nGFS 4 \
    --GridsPerR "${GRIDS_PER_R}" --ZMIN "${ZMIN}" --ZMAX "${ZMAX}" --RMAX "${RMAX}" \
    --folderToSave "${PILOT_ONE}/vel" --skip-video-encode
python3 "${WORK_DIR}/Video.py" --caseToProcess "${CASE_REAL}" --CPUs 4 --nGFS 4 \
    --GridsPerR "${GRIDS_PER_R}" --ZMIN "${ZMIN}" --ZMAX "${ZMAX}" --RMAX "${RMAX}" \
    --folderToSave "${PILOT_FOUR}/vel" --skip-video-encode
python3 "${WORK_DIR}/render_pair_video.py" "${CASE_REAL}" 10 /dev/null \
    "${ZMIN}" "${ZMAX}" "${RMAX}" "${FPS}" --cpus 1 --max-frames 4 \
    --frames-dir "${PILOT_ONE}/pair" --helper-dir "${WORK_DIR}" --skip-video-encode
python3 "${WORK_DIR}/render_pair_video.py" "${CASE_REAL}" 10 /dev/null \
    "${ZMIN}" "${ZMAX}" "${RMAX}" "${FPS}" --cpus 4 --max-frames 4 \
    --frames-dir "${PILOT_FOUR}/pair" --helper-dir "${WORK_DIR}" --skip-video-encode
diff -qr "${PILOT_ONE}/vel" "${PILOT_FOUR}/vel"
diff -qr "${PILOT_ONE}/pair" "${PILOT_FOUR}/pair"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
srun --exclusive --cpu-bind=cores --ntasks=1 --cpus-per-task="${THREADS}" \
    python3 "${WORK_DIR}/Video.py" --caseToProcess "${CASE_REAL}" --CPUs "${THREADS}" \
    --GridsPerR "${GRIDS_PER_R}" --ZMIN "${ZMIN}" --ZMAX "${ZMAX}" --RMAX "${RMAX}" \
    --folderToSave "${VEL_FRAMES}" --skip-video-encode
srun --exclusive --cpu-bind=cores --ntasks=1 --cpus-per-task="${THREADS}" \
    python3 "${WORK_DIR}/render_pair_video.py" "${CASE_REAL}" 10 /dev/null \
    "${ZMIN}" "${ZMAX}" "${RMAX}" "${FPS}" --cpus "${THREADS}" \
    --frames-dir "${PAIR_FRAMES}" --helper-dir "${WORK_DIR}" --skip-video-encode

VEL_COUNT="$(find "${VEL_FRAMES}" -maxdepth 1 -type f -name '*.png' -size +0c -printf . | wc -c)"
PAIR_COUNT="$(find "${PAIR_FRAMES}" -maxdepth 1 -type f -name '*.png' -size +0c -printf . | wc -c)"
[[ "${VEL_COUNT}" -eq "${SNAPSHOT_COUNT}" ]] || { echo "Velocity frame count mismatch" >&2; exit 4; }
[[ "${PAIR_COUNT}" -eq "${SNAPSHOT_COUNT}" ]] || { echo "Pair frame count mismatch" >&2; exit 4; }

printf 'case=%s\nsnapshots=%s\nvelocity_frames=%s\npair_frames=%s\ncompleted_at=%s\n' \
    "${CASE}" "${SNAPSHOT_COUNT}" "${VEL_COUNT}" "${PAIR_COUNT}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    > "${VOUT_REAL}/${CASE}_${TAG}_frames-complete.txt"
echo "FRAMES_COMPLETE case=${CASE} snapshots=${SNAPSHOT_COUNT}"
