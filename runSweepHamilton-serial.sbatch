#!/bin/bash
#SBATCH --job-name=burstBubble-S1
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH -p multi
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@durham.ac.uk
#SBATCH --output=slurm-stage1-%j.out
#SBATCH --error=slurm-stage1-%j.err

# ============================================================
# Bursting Bubble - Stage 1 Serial Runner (Hamilton)
# ============================================================
# This script runs Stage 1 (initialization) for all cases in parallel.
# Each case runs on 1 core (serial, no OpenMP/MPI).
# Up to 48 cases run simultaneously.
#
# Stage 1 generates restart files needed for Stage 2 (MPI simulation).
#
# AFTER RUNNING:
# Submit Stage 2 with: sbatch runSweepHamilton.sbatch
#
# SBATCH Parameters (Hamilton-specific):
#   -n 48: Number of tasks (one per case, max 48 concurrent)
#   --mem=64G: Memory (Stage 1 is less memory-intensive)
#   --time=12:00:00: Wall time (Stage 1 is quick per case)
#   -p multi: Partition for whole-node jobs
# ============================================================

set -euo pipefail  # Exit on error, unset variables, pipeline failures

# ============================================================
# Configuration
# ============================================================
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
SWEEP_FILE="${SCRIPT_DIR}/sweep.params"

# Maximum concurrent jobs (matches SLURM allocation)
MAX_CONCURRENT=${SLURM_NTASKS:-48}

# Stage 1 tmax - short run to generate restart file
STAGE1_TMAX="5e-2"

# ============================================================
# Validate Working Directory
# ============================================================
if [[ ! "$SCRIPT_DIR" =~ ^/nobackup/ ]]; then
    echo "WARNING: Running from $SCRIPT_DIR"
    echo "         Consider using /nobackup/$USER/ for better performance"
    echo "         Home directory has only 10GB quota and is slower"
    echo ""
fi

# ============================================================
# Print Job Information
# ============================================================
echo "============================================="
echo "Bursting Bubble - Stage 1 (Serial)"
echo "============================================="
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Max concurrent cases: ${MAX_CONCURRENT}"
echo "Partition: ${SLURM_JOB_PARTITION}"
echo ""
echo "Mode: Stage 1 only (serial initialization)"
echo "      Generates restart files for Stage 2"
echo "      Each case runs on 1 core (no OpenMP/MPI)"
echo ""

# ============================================================
# Load Required Modules (minimal for serial)
# ============================================================
echo "Loading modules..."
module purge
# No MPI needed for serial execution
echo "Modules loaded (serial mode - no MPI)"
echo ""

# ============================================================
# Setup Basilisk Environment
# ============================================================
echo "Setting up Basilisk environment..."

# Source version configuration
if [ -f "${SCRIPT_DIR}/src-local/basilisk_version.sh" ]; then
    # shellcheck disable=SC1091
    source "${SCRIPT_DIR}/src-local/basilisk_version.sh"
else
    echo "ERROR: src-local/basilisk_version.sh not found" >&2
    exit 1
fi

# Install Basilisk using centralized version (force reinstall)
install_basilisk --hard || exit 1
if [ -f "${SCRIPT_DIR}/.project_config" ]; then
    # shellcheck disable=SC1090
    source "${SCRIPT_DIR}/.project_config"
    echo "Basilisk environment loaded from .project_config"
    echo "BASILISK: $BASILISK"
else
    echo "ERROR: .project_config not found after installation" >&2
    exit 1
fi
echo ""

# ============================================================
# Validate Environment
# ============================================================
if [ -f "${SCRIPT_DIR}/src-local/parse_params.sh" ]; then
    # shellcheck disable=SC1091
    source "${SCRIPT_DIR}/src-local/parse_params.sh"
else
    echo "ERROR: src-local/parse_params.sh not found" >&2
    exit 1
fi

# Source sweep utilities library
if [ -f "${SCRIPT_DIR}/src-local/sweep_utils.sh" ]; then
    # shellcheck disable=SC1091
    source "${SCRIPT_DIR}/src-local/sweep_utils.sh"
else
    echo "ERROR: src-local/sweep_utils.sh not found" >&2
    exit 1
fi

if [ ! -f "$SWEEP_FILE" ]; then
    echo "ERROR: Sweep file not found: $SWEEP_FILE" >&2
    exit 1
fi

echo "Sweep file: $SWEEP_FILE"
echo ""

# ============================================================
# Parse Sweep Configuration
# ============================================================
echo "Parsing sweep configuration..."

# shellcheck disable=SC1090
source "$SWEEP_FILE"

if [ -z "${BASE_CONFIG:-}" ]; then
    echo "ERROR: BASE_CONFIG not defined in sweep file" >&2
    exit 1
fi

if [ ! -f "$BASE_CONFIG" ]; then
    echo "ERROR: Base configuration file not found: $BASE_CONFIG" >&2
    exit 1
fi

# Validate case range using shared utility
validate_case_range "$CASE_START" "$CASE_END" || exit 1

echo "Base configuration: $BASE_CONFIG"
echo "Case number range: $CASE_START to $CASE_END"
echo ""

# ============================================================
# Extract Sweep Variables (using shared utility)
# ============================================================
extract_sweep_variables "$SWEEP_FILE" || exit 1
print_sweep_variables
echo ""

# ============================================================
# Generate Parameter Combinations (using shared utility)
# ============================================================
echo "Generating parameter combinations..."

# Use work directory for temp files (HPC: avoid slow /tmp)
setup_sweep_temp_dir "$SCRIPT_DIR" || exit 1

# Generate combinations with verbose output
generate_sweep_combinations 1 || exit 1

echo "Generated $SWEEP_COMBINATION_COUNT parameter combinations"

# Validate combination count
validate_combination_count "$CASE_START" "$CASE_END" || exit 1

echo ""

# ============================================================
# Prepare Status Tracking
# ============================================================
STATUS_DIR="${SWEEP_TEMP_DIR}/status"
mkdir -p "$STATUS_DIR"

# ============================================================
# Function to Run Single Case (Stage 1)
# ============================================================
run_stage1_case() {
    local param_file=$1
    local status_dir=$2
    local script_dir=$3
    local stage1_tmax=$4

    # Source parameter parsing in subshell
    source "${script_dir}/src-local/parse_params.sh"

    parse_param_file "$param_file"
    local CASE_NO=$(get_param "CaseNo")
    local Oh=$(get_param "Oh" "1e-2")
    local Bond=$(get_param "Bond" "1e-3")
    local MAXlevel=$(get_param "MAXlevel" "10")
    local zWall=$(get_param "zWall" "4")

    if [ -z "$CASE_NO" ]; then
        echo "FAILED" > "${status_dir}/${CASE_NO:-unknown}"
        return 1
    fi

    local CASE_DIR="${script_dir}/simulationCases/${CASE_NO}"
    local LOG_FILE="${CASE_DIR}/stage1.log"

    # Create case directory BEFORE redirecting to log file
    mkdir -p "$CASE_DIR"

    # Use existing case.params if present (allows manual parameter edits for reruns)
    # Otherwise, create from sweep-generated params
    if [ ! -f "$CASE_DIR/case.params" ]; then
        cp "$param_file" "$CASE_DIR/case.params"
    fi
    # Re-parse from case.params (source of truth for this case)
    parse_param_file "$CASE_DIR/case.params"
    Oh=$(get_param "Oh" "1e-2")
    Bond=$(get_param "Bond" "1e-3")
    MAXlevel=$(get_param "MAXlevel" "10")
    zWall=$(get_param "zWall" "4")

    {
        echo "========================================="
        echo "Case $CASE_NO - Stage 1 ($(date))"
        echo "========================================="
        echo "Parameters: Oh=$Oh, Bond=$Bond"
        echo "            MAXlevel=$MAXlevel, tmax=$stage1_tmax, zWall=$zWall"

        cd "$CASE_DIR"

        # Copy source file
        local SRC_FILE_ORIG="${script_dir}/simulationCases/burstingBubble.c"
        if [ ! -f "$SRC_FILE_ORIG" ]; then
            echo "ERROR: Source file not found"
            echo "FAILED" > "${status_dir}/${CASE_NO}"
            return 1
        fi
        # Use existing source file if present (allows local code edits for reruns)
        if [ ! -f "burstingBubble.c" ]; then
            cp "$SRC_FILE_ORIG" "burstingBubble.c"
        fi

        # Create symlink to DataFiles
        if [ ! -e "DataFiles" ]; then
            ln -s ../DataFiles DataFiles
        fi

        # Compile SERIAL (no OpenMP, no MPI)
        echo "Compiling (serial)..."
        if ! qcc -I../../src-local -Wall -O2 -disable-dimensions \
            burstingBubble.c -o burstingBubble -lm 2>&1; then
            echo "ERROR: Compilation failed"
            echo "FAILED" > "${status_dir}/${CASE_NO}"
            return 1
        fi
        echo "Compilation successful"

        # Run Stage 1 with small tmax
        # Bursting-Bubble argument order: $MAXlevel $Oh $Bond $tmax $zWall
        echo "Running Stage 1 (tmax=$stage1_tmax)..."
        if ./burstingBubble $MAXlevel $Oh $Bond $stage1_tmax $zWall; then
            echo "Stage 1 completed"

            # Verify restart file was created
            if [ -f "restart" ]; then
                echo "Restart file created successfully"
                echo "SUCCESS" > "${status_dir}/${CASE_NO}"
            else
                echo "WARNING: restart file not found after Stage 1"
                echo "FAILED" > "${status_dir}/${CASE_NO}"
                return 1
            fi
        else
            echo "ERROR: Stage 1 failed"
            echo "FAILED" > "${status_dir}/${CASE_NO}"
            return 1
        fi

        echo "Case $CASE_NO Stage 1 complete"
    } > "$LOG_FILE" 2>&1

    return 0
}

export -f run_stage1_case

# ============================================================
# Run All Cases in Parallel (Stage 1)
# ============================================================
echo "============================================="
echo "Running $SWEEP_COMBINATION_COUNT Cases (Stage 1)"
echo "============================================="
echo "Max concurrent: ${MAX_CONCURRENT}"
echo "Each case: serial execution, tmax=${STAGE1_TMAX}"
echo ""

RUNNING_JOBS=0
PIDS=()

for param_file in "${SWEEP_CASE_FILES[@]}"; do
    # Wait if we've hit the concurrency limit
    while [ $RUNNING_JOBS -ge $MAX_CONCURRENT ]; do
        # Wait for any background job to finish
        wait -n 2>/dev/null || true
        ((RUNNING_JOBS--)) || true
    done

    # Launch case in background
    run_stage1_case "$param_file" "$STATUS_DIR" "$SCRIPT_DIR" "$STAGE1_TMAX" &
    PIDS+=($!)
    ((RUNNING_JOBS++)) || true

    # Extract case number for logging
    parse_param_file "$param_file"
    CASE_NO=$(get_param "CaseNo")
    echo "Launched case $CASE_NO (PID: ${PIDS[-1]}, running: $RUNNING_JOBS)"
done

echo ""
echo "All cases launched, waiting for completion..."
echo ""

# Wait for all remaining jobs
wait || true

# ============================================================
# Collect Results
# ============================================================
SUCCESSFUL_CASES=0
FAILED_CASES=0

for status_file in "${STATUS_DIR}"/*; do
    if [ -f "$status_file" ]; then
        status=$(cat "$status_file")
        case_no=$(basename "$status_file")
        if [ "$status" = "SUCCESS" ]; then
            ((SUCCESSFUL_CASES++)) || true
        else
            ((FAILED_CASES++)) || true
            echo "FAILED: Case $case_no"
        fi
    fi
done

# ============================================================
# Final Summary
# ============================================================
echo ""
echo "============================================="
echo "Stage 1 Complete"
echo "============================================="
echo "Job completed at: $(date)"
echo "Total cases: $SWEEP_COMBINATION_COUNT"
echo "Successful: $SUCCESSFUL_CASES"
echo "Failed: $FAILED_CASES"
echo ""
echo "Restart files location: simulationCases/*/restart"
echo ""
if [ $SUCCESSFUL_CASES -gt 0 ]; then
    echo "Next step: Submit Stage 2 with:"
    echo "  sbatch runSweepHamilton.sbatch"
fi
echo "============================================="

if [ $FAILED_CASES -gt 0 ]; then
    echo "WARNING: $FAILED_CASES case(s) failed. Check logs in simulationCases/*/stage1.log" >&2
    exit 1
fi

exit 0
