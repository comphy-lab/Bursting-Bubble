#!/bin/bash
#SBATCH --job-name=burstingBubbleSweep
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --time=100:00:00
#SBATCH --partition=genoa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vatsal.sanjay@comphy-lab.org
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# ============================================================
# Bursting Bubble - HPC Parameter Sweep Runner (Snellius)
# ============================================================
# This script runs Stage 2 (MPI simulation) on Snellius HPC.
#
# PREREQUISITE: Restart files must already exist for each case.
# Run Stage 1 locally before submitting:
#   ./runParameterSweep.sh --stage1-only sweep.params
#
# BEFORE RUNNING:
# 1. Run Stage 1 locally to generate restart files for all cases
# 2. Install Basilisk locally (creates basilisk/ and .project_config):
#    curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-01-13
# 3. Edit sweep.params to set parameter ranges and case numbers
# 4. Adjust SBATCH parameters above (especially --ntasks, --time)
# 5. Submit with: sbatch runSweepSnellius.sbatch
#
# SBATCH Parameters to Customize:
#   --ntasks: Number of MPI tasks per case (currently 48)
#   --time: Wall time for entire sweep (currently 100 hours)
#   --partition: Compute partition (currently genoa)
#   --mail-user: Email for job notifications
# ============================================================

set -euo pipefail

# ============================================================
# Configuration
# ============================================================
# Use SLURM_SUBMIT_DIR (directory where sbatch was called) as the project root
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)"
SWEEP_FILE="${SCRIPT_DIR}/sweep.params"

# Source file to compile inside each case directory.
# Update for other cases (e.g., asymmetric coalescence) if needed.
SOURCE_FILE_NAME="burstingBubble.c"
EXECUTABLE_NAME="${SOURCE_FILE_NAME%.c}"

# ============================================================
# Print Job Information
# ============================================================
echo "============================================="
echo "Bursting Bubble - HPC Parameter Sweep"
echo "============================================="
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Number of MPI tasks: ${SLURM_NTASKS:-unknown}"
echo "Partition: ${SLURM_JOB_PARTITION:-unknown}"
echo ""
echo "Mode: Stage 2 only (MPI simulation)"
echo "      Restart files must exist for each case"
echo ""

# ============================================================
# Load Required Modules
# ============================================================
echo "Loading modules..."
module purge
module load 2024
module load OpenMPI/5.0.3-GCC-13.3.0
echo "Modules loaded successfully"
echo ""

# ============================================================
# Setup Basilisk Environment
# ============================================================
echo "Setting up Basilisk environment..."
cd "$SCRIPT_DIR"

# Source version configuration
if [ -f "${SCRIPT_DIR}/src-local/basilisk_version.sh" ]; then
  # shellcheck disable=SC1091
  source "${SCRIPT_DIR}/src-local/basilisk_version.sh"
else
  echo "ERROR: src-local/basilisk_version.sh not found" >&2
  exit 1
fi

# Install Basilisk using centralized version
install_basilisk || exit 1

if [ -f "${SCRIPT_DIR}/.project_config" ]; then
  # shellcheck disable=SC1090
  source "${SCRIPT_DIR}/.project_config"
  echo "Basilisk environment loaded from .project_config"
  echo "BASILISK: ${BASILISK:-unset}"
else
  echo "ERROR: .project_config not found after Basilisk install" >&2
  exit 1
fi
echo ""

# ============================================================
# Validate Environment
# ============================================================
# Source parameter parsing library
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
# shellcheck disable=SC1090
source "$SWEEP_FILE"

# Validate required variables
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
# Run Simulations (Stage 2 Only - MPI)
# ============================================================
echo "============================================="
echo "Running $SWEEP_COMBINATION_COUNT Simulations"
echo "============================================="
echo "Stage 2 only: MPI simulation (${SLURM_NTASKS:-unknown} tasks)"
echo "Requires existing restart files in each case directory"
echo ""

# Counters for tracking
SUCCESSFUL_CASES=0
FAILED_CASES=0
SKIPPED_CASES=0

# Initialize progress tracking
sweep_progress_init "$SWEEP_COMBINATION_COUNT"

for param_file in "${SWEEP_CASE_FILES[@]}"; do
  # Parse parameter file to get values
  parse_param_file "$param_file"
  CASE_NO=$(get_param "CaseNo")
  Oh=$(get_param "Oh" "1e-2")
  Bond=$(get_param "Bond" "1e-3")
  MAXlevel=$(get_param "MAXlevel" "10")
  tmax=$(get_param "tmax" "1.0")
  zWall=$(get_param "zWall" "0.025")

  if [ -z "${CASE_NO:-}" ]; then
    echo "ERROR: CaseNo not found in $param_file" >&2
    FAILED_CASES=$((FAILED_CASES + 1))
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
    SKIPPED_CASES=$((SKIPPED_CASES + 1))
    continue
  fi

  # Change to case directory
  cd "$CASE_DIR"

  # Check for restart file (CRITICAL)
  if [ ! -f "restart" ]; then
    echo "ERROR: restart file not found for case $CASE_NO" >&2
    echo "       Run Stage 1 locally first:" >&2
    echo "       ./runSimulation.sh --stage1 <params_file>" >&2
    SKIPPED_CASES=$((SKIPPED_CASES + 1))
    cd "$SCRIPT_DIR"
    continue
  fi

  echo "Found restart file"

  # Use existing case.params if present (allows manual parameter edits for reruns)
  if [ -f "$CASE_DIR/case.params" ]; then
    echo "Using existing case.params (manual edits preserved)"
    parse_param_file "$CASE_DIR/case.params"
    Oh=$(get_param "Oh" "1e-2")
    Bond=$(get_param "Bond" "1e-3")
    MAXlevel=$(get_param "MAXlevel" "10")
    tmax=$(get_param "tmax" "1.0")
    zWall=$(get_param "zWall" "0.025")
  fi

  # Copy source file to case directory if needed
  SRC_FILE_ORIG="${SCRIPT_DIR}/simulationCases/${SOURCE_FILE_NAME}"
  SRC_FILE_LOCAL="$CASE_DIR/$SOURCE_FILE_NAME"

  if [ ! -f "$SRC_FILE_ORIG" ]; then
    echo "ERROR: Source file $SRC_FILE_ORIG not found" >&2
    FAILED_CASES=$((FAILED_CASES + 1))
    cd "$SCRIPT_DIR"
    continue
  fi

  # Use existing source file if present (allows local code edits for reruns)
  if [ ! -f "$SRC_FILE_LOCAL" ]; then
    cp "$SRC_FILE_ORIG" "$SRC_FILE_LOCAL"
    echo "Copied source file to case directory"
  else
    echo "Using existing $SOURCE_FILE_NAME (local edits preserved)"
  fi

  # Create symlink to DataFiles if needed
  if [ ! -e "DataFiles" ]; then
    ln -s ../DataFiles DataFiles
  fi

  # ============================================================
  # Stage 2: Full Simulation with MPI
  # ============================================================
  echo ""
  echo "Running Stage 2: Full simulation (MPI)..."

  # Compile with MPI
  if CC99='mpicc -std=c99 -D_GNU_SOURCE=1' qcc \
    -Wall -O2 -D_MPI=1 -disable-dimensions \
    "$SOURCE_FILE_NAME" -o "$EXECUTABLE_NAME" -lm 2>&1; then
    echo "MPI compilation successful"
  else
    echo "ERROR: MPI compilation failed for case $CASE_NO" >&2
    FAILED_CASES=$((FAILED_CASES + 1))
    cd "$SCRIPT_DIR"
    continue
  fi

  # Run simulation with srun
  echo "Running full simulation with ${SLURM_NTASKS:-unknown} MPI tasks..."
  echo "Command: srun -n ${SLURM_NTASKS:-1} ./$EXECUTABLE_NAME $MAXlevel $Oh $Bond $tmax $zWall"
  echo ""

  if srun -n "${SLURM_NTASKS:-1}" ./$EXECUTABLE_NAME $MAXlevel $Oh $Bond $tmax $zWall; then
    echo ""
    echo "Case $CASE_NO completed successfully"
    SUCCESSFUL_CASES=$((SUCCESSFUL_CASES + 1))
  else
    EXIT_CODE=$?
    echo ""
    echo "ERROR: Case $CASE_NO failed with exit code $EXIT_CODE" >&2
    FAILED_CASES=$((FAILED_CASES + 1))
  fi

  # Return to root directory
  cd "$SCRIPT_DIR"

  # Update progress
  sweep_progress_update "$CASE_NO"
  echo ""
done

# ============================================================
# Final Summary
# ============================================================
echo "============================================="
echo "Parameter Sweep Complete"
echo "============================================="
echo "Job completed at: $(date)"
echo "Total cases: $SWEEP_COMBINATION_COUNT"
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
