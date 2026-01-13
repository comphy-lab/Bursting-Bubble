#!/bin/bash
#SBATCH --job-name=burstingBubble
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -c 1
#SBATCH --mem=240G
#SBATCH --time=3-00:00:00
#SBATCH --gres=tmp:100G
#SBATCH -p multi
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@durham.ac.uk
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# ============================================================
# Bursting Bubble - HPC Parameter Sweep Runner (Hamilton)
# ============================================================
# This script runs Stage 2 (MPI simulation) on Durham University's Hamilton HPC.
#
# PREREQUISITE: Restart files must already exist for each case.
# Run Stage 1 locally before submitting:
#   ./runParameterSweep.sh --stage1-only sweep.params
#
# BEFORE RUNNING:
# 1. Run Stage 1 locally to generate restart files for all cases
# 2. Clone repo to /nobackup/$USER/ (NOT home directory)
# 3. Install Basilisk locally (creates basilisk/ and .project_config):
#    curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-01-13
# 4. Edit sweep.params to set parameter ranges and case numbers
# 5. Update --mail-user above with your Durham email
# 6. Submit with: sbatch runSweepHamilton.sbatch
#
# SBATCH Parameters (Hamilton-specific):
#   -n 128: Number of MPI tasks (whole node = 128 cores)
#   --mem=240G: Memory per node (max 250G on standard nodes)
#   --time=3-00:00:00: Max wall time (3 days for multi partition)
#   --gres=tmp:100G: Temporary disk space ($TMPDIR)
#   -p multi: Partition for whole-node jobs
# ============================================================

set -e  # Exit on error (but we'll handle case failures gracefully)
set -o pipefail

# ============================================================
# Configuration
# ============================================================
# Use SLURM_SUBMIT_DIR (directory where sbatch was called) as the project root
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
SWEEP_FILE="${SCRIPT_DIR}/sweep.params"

# Source file to compile inside each case directory
SOURCE_FILE_NAME="burstingBubble.c"
EXECUTABLE_NAME="${SOURCE_FILE_NAME%.c}"

# Number of MPI tasks for Stage 2
MPI_TASKS=${SLURM_NTASKS:-128}

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
echo "Bursting Bubble - HPC Parameter Sweep"
echo "============================================="
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Number of MPI tasks: ${MPI_TASKS}"
echo "Partition: ${SLURM_JOB_PARTITION}"
echo "Memory allocated: ${SLURM_MEM_PER_NODE}MB"
echo ""
echo "Mode: Stage 2 only (MPI simulation)"
echo "      Restart files must exist for each case"
echo ""

# ============================================================
# Load Required Modules
# ============================================================
echo "Loading modules..."
module purge
module load gcc
module load openmpi
echo "Modules loaded successfully"
echo ""

# ============================================================
# Setup Basilisk Environment
# ============================================================
echo "Setting up Basilisk environment..."
echo "Ensuring Basilisk is installed locally (ref v2026-01-13)..."
cd "$SCRIPT_DIR"

if ! command -v curl >/dev/null 2>&1; then
    echo "ERROR: curl not found (required to install Basilisk)" >&2
    exit 1
fi

curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-01-13

if [ -f "${SCRIPT_DIR}/.project_config" ]; then
    source "${SCRIPT_DIR}/.project_config"
    echo "Basilisk environment loaded from .project_config"
    echo "BASILISK: $BASILISK"
else
    echo "ERROR: .project_config not found after Basilisk install" >&2
    echo "       Basilisk install likely failed; see slurm output for details" >&2
    exit 1
fi
echo ""

# ============================================================
# Validate Environment
# ============================================================
# Source parameter parsing library
if [ -f "${SCRIPT_DIR}/src-local/parse_params.sh" ]; then
    source "${SCRIPT_DIR}/src-local/parse_params.sh"
else
    echo "ERROR: src-local/parse_params.sh not found" >&2
    exit 1
fi

# Check sweep file exists
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

# Source the sweep file to get variables
source "$SWEEP_FILE"

# Validate required variables
if [ -z "$BASE_CONFIG" ]; then
    echo "ERROR: BASE_CONFIG not defined in sweep file" >&2
    exit 1
fi

if [ -z "$CASE_START" ] || [ -z "$CASE_END" ]; then
    echo "ERROR: CASE_START and CASE_END must be defined in sweep file" >&2
    exit 1
fi

# Validate CaseNo range
if [ "$CASE_START" -lt 1000 ] || [ "$CASE_START" -gt 9999 ]; then
    echo "ERROR: CASE_START must be 4-digit (1000-9999), got: $CASE_START" >&2
    exit 1
fi

if [ "$CASE_END" -lt "$CASE_START" ] || [ "$CASE_END" -gt 9999 ]; then
    echo "ERROR: CASE_END must be >= CASE_START and <= 9999, got: $CASE_END" >&2
    exit 1
fi

if [ ! -f "$BASE_CONFIG" ]; then
    echo "ERROR: Base configuration file not found: $BASE_CONFIG" >&2
    exit 1
fi

echo "Base configuration: $BASE_CONFIG"
echo "Case number range: $CASE_START to $CASE_END"
echo ""

# ============================================================
# Extract Sweep Variables
# ============================================================
SWEEP_VARS=()
SWEEP_VALUES=()

# Read sweep file and extract SWEEP_* variables
while IFS='=' read -r key value; do
    # Skip comments and empty lines
    [[ "$key" =~ ^[[:space:]]*# ]] && continue
    [[ -z "$key" ]] && continue

    # Match SWEEP_* variables
    if [[ "$key" =~ ^[[:space:]]*SWEEP_([^=]+) ]]; then
        var_name="${BASH_REMATCH[1]}"
        # Remove inline comments and whitespace
        value=$(echo "$value" | sed 's/#.*//' | xargs)

        SWEEP_VARS+=("$var_name")
        SWEEP_VALUES+=("$value")
    fi
done < "$SWEEP_FILE"

if [ ${#SWEEP_VARS[@]} -eq 0 ]; then
    echo "ERROR: No SWEEP_* variables found in $SWEEP_FILE" >&2
    exit 1
fi

echo "Sweep variables:"
for i in "${!SWEEP_VARS[@]}"; do
    echo "  ${SWEEP_VARS[$i]} = ${SWEEP_VALUES[$i]}"
done
echo ""

# ============================================================
# Generate Parameter Combinations
# ============================================================
echo "Generating parameter combinations..."

# Create temporary directory for generated parameter files
# Use work directory, not /tmp (node-local on Hamilton)
TEMP_DIR="${SCRIPT_DIR}/.sweep_tmp_$$"
mkdir -p "$TEMP_DIR" || {
    echo "ERROR: Failed to create temp directory: $TEMP_DIR" >&2
    exit 1
}
trap "rm -rf $TEMP_DIR" EXIT

CASE_NUM=$CASE_START
COMBINATION_COUNT=0
CASE_FILES=()

# Recursive function to generate all combinations
generate_combinations() {
    local depth=$1
    shift
    local current_values=("$@")

    if [ $depth -eq ${#SWEEP_VARS[@]} ]; then
        # Base case: all variables assigned, create parameter file
        local case_file="${TEMP_DIR}/case_$(printf "%04d" $CASE_NUM).params"

        # Copy base config
        cp "$BASE_CONFIG" "$case_file"

        # Override CaseNo
        if grep -q "^CaseNo=" "$case_file"; then
            sed -i'.bak' "s|^CaseNo=.*|CaseNo=${CASE_NUM}|" "$case_file"
        else
            echo "CaseNo=${CASE_NUM}" >> "$case_file"
        fi
        rm -f "${case_file}.bak"

        # Override with sweep values
        for i in "${!SWEEP_VARS[@]}"; do
            local var="${SWEEP_VARS[$i]}"
            local val="${current_values[$i]}"

            if grep -q "^${var}=" "$case_file"; then
                sed -i'.bak' "s|^${var}=.*|${var}=${val}|" "$case_file"
            else
                echo "${var}=${val}" >> "$case_file"
            fi
            rm -f "${case_file}.bak"
        done

        CASE_FILES+=("$case_file")
        ((COMBINATION_COUNT++)) || true

        # Print summary
        echo "Case $CASE_NUM:"
        for i in "${!SWEEP_VARS[@]}"; do
            echo "  ${SWEEP_VARS[$i]} = ${current_values[$i]}"
        done
        echo ""

        ((CASE_NUM++)) || true
        return
    fi

    # Recursive case: iterate through values for current variable
    local values="${SWEEP_VALUES[$depth]}"
    IFS=',' read -ra value_array <<< "$values"

    for val in "${value_array[@]}"; do
        val=$(echo "$val" | xargs)  # Trim whitespace
        generate_combinations $((depth + 1)) "${current_values[@]}" "$val"
    done
}

# Start recursion
generate_combinations 0

echo "Generated $COMBINATION_COUNT parameter combinations"

# Check if number of combinations matches the range
EXPECTED_COUNT=$((CASE_END - CASE_START + 1))
if [ $COMBINATION_COUNT -ne $EXPECTED_COUNT ]; then
    echo "WARNING: Generated $COMBINATION_COUNT combinations, but CASE_END suggests $EXPECTED_COUNT" >&2
    echo "         Consider adjusting CASE_END in sweep file" >&2
fi

if [ $COMBINATION_COUNT -gt $EXPECTED_COUNT ]; then
    echo "ERROR: Too many combinations ($COMBINATION_COUNT) for range $CASE_START-$CASE_END" >&2
    exit 1
fi

echo ""

# ============================================================
# Run Simulations (Stage 2 Only - MPI)
# ============================================================
echo "============================================="
echo "Running $COMBINATION_COUNT Simulations"
echo "============================================="
echo "Stage 2 only: MPI simulation (${MPI_TASKS} tasks)"
echo "Requires existing restart files in each case directory"
echo ""

# Counters for tracking
SUCCESSFUL_CASES=0
FAILED_CASES=0
SKIPPED_CASES=0

for param_file in "${CASE_FILES[@]}"; do
    # Parse parameter file to get values
    parse_param_file "$param_file"
    CASE_NO=$(get_param "CaseNo")
    MAXlevel=$(get_param "MAXlevel" "10")
    Oh=$(get_param "Oh" "1e-2")
    Bond=$(get_param "Bond" "1e-3")
    tmax=$(get_param "tmax" "1.0")
    zWall=$(get_param "zWall" "0.025")

    if [ -z "$CASE_NO" ]; then
        echo "ERROR: CaseNo not found in $param_file" >&2
        ((FAILED_CASES++)) || true
        continue
    fi

    CASE_DIR="${SCRIPT_DIR}/simulationCases/${CASE_NO}"

    echo "========================================="
    echo "Case $CASE_NO ($(date))"
    echo "========================================="
    echo "Parameters: Oh=$Oh, Bond=$Bond"
    echo "            MAXlevel=$MAXlevel, tmax=$tmax, zWall=$zWall"

    # Check if case directory exists
    if [ ! -d "$CASE_DIR" ]; then
        echo "ERROR: Case directory not found: $CASE_DIR" >&2
        echo "       Run Stage 1 locally first" >&2
        ((SKIPPED_CASES++)) || true
        continue
    fi

    # Change to case directory
    cd "$CASE_DIR"

    # Check for restart file (CRITICAL)
    if [ ! -f "restart" ]; then
        echo "ERROR: restart file not found for case $CASE_NO" >&2
        echo "       Run Stage 1 locally first:" >&2
        echo "       ./runSimulation.sh --stage1 <params_file>" >&2
        ((SKIPPED_CASES++)) || true
        cd "$SCRIPT_DIR"
        continue
    fi

    echo "Found restart file"

    # Copy source file to case directory if needed
    SRC_FILE_ORIG="${SCRIPT_DIR}/simulationCases/${SOURCE_FILE_NAME}"
    SRC_FILE_LOCAL="$CASE_DIR/$SOURCE_FILE_NAME"

    if [ ! -f "$SRC_FILE_ORIG" ]; then
        echo "ERROR: Source file $SRC_FILE_ORIG not found" >&2
        ((FAILED_CASES++)) || true
        cd "$SCRIPT_DIR"
        continue
    fi

    cp "$SRC_FILE_ORIG" "$SRC_FILE_LOCAL"

    # Create symlink to DataFiles if needed
    if [ ! -e "DataFiles" ]; then
        ln -s ../DataFiles DataFiles
    fi

    # ============================================================
    # Stage 2: Full Simulation with MPI
    # ============================================================
    echo ""
    echo "Stage 2: Full simulation (MPI)..."

    # Compile with MPI
    if CC99='mpicc -std=c99 -D_GNU_SOURCE=1' qcc -I../../src-local \
        -Wall -O2 -D_MPI=1 -disable-dimensions \
        "$SOURCE_FILE_NAME" -o "$EXECUTABLE_NAME" -lm 2>&1; then
        echo "MPI compilation successful"
    else
        echo "ERROR: MPI compilation failed for case $CASE_NO" >&2
        ((FAILED_CASES++)) || true
        cd "$SCRIPT_DIR"
        continue
    fi

    # Run simulation with mpirun (Hamilton uses mpirun, not srun)
    echo "Running full simulation with ${MPI_TASKS} MPI tasks..."
    echo "Command: mpirun -np ${MPI_TASKS} ./$EXECUTABLE_NAME $MAXlevel $Oh $Bond $tmax $zWall"
    echo ""

    if mpirun -np ${MPI_TASKS} ./$EXECUTABLE_NAME $MAXlevel $Oh $Bond $tmax $zWall; then
        echo ""
        echo "Case $CASE_NO completed successfully"
        ((SUCCESSFUL_CASES++)) || true
    else
        EXIT_CODE=$?
        echo ""
        echo "ERROR: Case $CASE_NO failed with exit code $EXIT_CODE" >&2
        ((FAILED_CASES++)) || true
    fi

    # Return to root directory
    cd "$SCRIPT_DIR"
    echo ""
done

# ============================================================
# Final Summary
# ============================================================
echo "============================================="
echo "Parameter Sweep Complete"
echo "============================================="
echo "Job completed at: $(date)"
echo "Total cases: $COMBINATION_COUNT"
echo "Successful: $SUCCESSFUL_CASES"
echo "Failed: $FAILED_CASES"
echo "Skipped (no restart file): $SKIPPED_CASES"
echo "Output location: simulationCases/"
echo "============================================="
echo ""

# Exit with error if any cases failed or were skipped
if [ $FAILED_CASES -gt 0 ] || [ $SKIPPED_CASES -gt 0 ]; then
    echo "WARNING: $FAILED_CASES case(s) failed, $SKIPPED_CASES case(s) skipped." >&2
    echo "         Check logs for details." >&2
    exit 1
fi

exit 0
