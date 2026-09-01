#!/bin/bash
#SBATCH --job-name=sbbvidencode
#SBATCH --account=nctth181
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --partition=genoa
#SBATCH --output=slurm-video-encode-%x-%j.out
#SBATCH --error=slurm-video-encode-%x-%j.err
# =============================================================================
# Encode and verify the two frame products emitted by runRenderVideos.sbatch.
# Submit with --dependency=afterok:<frame-job-id> so a small allocation, rather
# than a 192-core node, owns serial assembly and checksums.
#
# Required environment: CASE, CASEDIR, TAG, VOUT.
# Optional: FRAMING="ZMIN ZMAX RMAX FPS", FIELD_INPUT_FPS (90),
# FIELD_OUTPUT_FPS (30).
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
VOUT_REAL="$(realpath "${VOUT}")"
[[ "${CASE_REAL}" == "${SUBMIT_REAL}"/simulationCases/"${CASE}" ]] || exit 2
[[ "${VOUT_REAL}" == "${SUBMIT_REAL}"/postprocessed/* ]] || exit 2

read -r _ _ _ PAIR_FPS <<< "${FRAMING:--2.05 1.5 1.1 18}"
FIELD_INPUT_FPS="${FIELD_INPUT_FPS:-90}"
FIELD_OUTPUT_FPS="${FIELD_OUTPUT_FPS:-30}"
for value in "${PAIR_FPS}" "${FIELD_INPUT_FPS}" "${FIELD_OUTPUT_FPS}"; do
    [[ "${value}" =~ ^[0-9]+$ ]] && [[ "${value}" -gt 0 ]] || exit 2
done

module purge
module load 2024
module load FFmpeg/7.0.2-GCCcore-13.3.0
command -v ffmpeg >/dev/null
command -v ffprobe >/dev/null

MARKER="${VOUT_REAL}/${CASE}_${TAG}_frames-complete.txt"
test -s "${MARKER}"
EXPECTED="$(awk -F= '$1=="snapshots" {print $2}' "${MARKER}")"
[[ "${EXPECTED}" =~ ^[0-9]+$ ]] && [[ "${EXPECTED}" -gt 0 ]] || exit 3
VEL_FRAMES="${VOUT_REAL}/frames_${CASE}_vel"
PAIR_FRAMES="${VOUT_REAL}/frames_${CASE}_pair"
VEL_COUNT="$(find "${VEL_FRAMES}" -maxdepth 1 -type f -name '*.png' -size +0c -printf . | wc -c)"
PAIR_COUNT="$(find "${PAIR_FRAMES}" -maxdepth 1 -type f -name '*.png' -size +0c -printf . | wc -c)"
[[ "${VEL_COUNT}" -eq "${EXPECTED}" ]] || { echo "Velocity frame count mismatch" >&2; exit 3; }
[[ "${PAIR_COUNT}" -eq "${EXPECTED}" ]] || { echo "Pair frame count mismatch" >&2; exit 3; }

VEL_OUT="${VOUT_REAL}/${CASE}_${TAG}_vel_dissipation.mp4"
PAIR_OUT="${VOUT_REAL}/${CASE}_${TAG}_interface_marker_AMR.mp4"
VEL_TMP="${VOUT_REAL}/.${CASE}_${TAG}_vel_dissipation.${SLURM_JOB_ID}.tmp.mp4"
PAIR_TMP="${VOUT_REAL}/.${CASE}_${TAG}_interface_marker_AMR.${SLURM_JOB_ID}.tmp.mp4"
trap 'rm -f -- "${VEL_TMP}" "${PAIR_TMP}"' EXIT

ffmpeg -y -threads "${SLURM_CPUS_PER_TASK:-4}" -framerate "${FIELD_INPUT_FPS}" \
    -pattern_type glob -i "${VEL_FRAMES}/*.png" -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' \
    -c:v libx264 -r "${FIELD_OUTPUT_FPS}" -pix_fmt yuv420p "${VEL_TMP}"
ffmpeg -y -threads "${SLURM_CPUS_PER_TASK:-4}" -framerate "${PAIR_FPS}" \
    -i "${PAIR_FRAMES}/%06d.png" -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' \
    -c:v libx264 -pix_fmt yuv420p -crf 20 "${PAIR_TMP}"

for video in "${VEL_TMP}" "${PAIR_TMP}"; do
    test -s "${video}"
    ffprobe -v error -select_streams v:0 -show_entries stream=codec_name,width,height,duration \
        -of default=noprint_wrappers=1 "${video}"
done
mv "${VEL_TMP}" "${VEL_OUT}"
mv "${PAIR_TMP}" "${PAIR_OUT}"
(
    cd "${VOUT_REAL}"
    sha256sum "$(basename "${VEL_OUT}")" "$(basename "${PAIR_OUT}")" \
        > "${CASE}_${TAG}_videos.sha256"
)
printf 'case=%s\nsnapshots=%s\nvelocity_video=%s\npair_video=%s\ncompleted_at=%s\n' \
    "${CASE}" "${EXPECTED}" "${VEL_OUT}" "${PAIR_OUT}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    > "${VOUT_REAL}/${CASE}_${TAG}_videos-complete.txt"
trap - EXIT
echo "VIDEOS_COMPLETE case=${CASE} snapshots=${EXPECTED}"
