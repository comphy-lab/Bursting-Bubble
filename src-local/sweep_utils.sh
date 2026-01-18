#!/bin/bash
# sweep_utils.sh - Shared utility functions for parameter sweep generation
#
# Description:
#   Provides reusable functions for extracting sweep parameters, validating
#   case ranges, and generating Cartesian product combinations of parameters.
#   This library is sourced by sweep scripts to eliminate code duplication.
#
# Functions:
#   extract_sweep_variables <file>     - Extract SWEEP_* variables from file
#   validate_case_range <start> <end>  - Validate 4-digit case number range
#   generate_sweep_combinations        - Generate all parameter combinations
#   setup_sweep_temp_dir [base_dir]    - Create temp directory with cleanup
#
# Usage:
#   source src-local/sweep_utils.sh
#   source "$SWEEP_FILE"  # Get BASE_CONFIG, CASE_START, CASE_END
#   extract_sweep_variables "$SWEEP_FILE"
#   validate_case_range "$CASE_START" "$CASE_END"
#   setup_sweep_temp_dir "$SCRIPT_DIR"
#   generate_sweep_combinations
#
# Dependencies:
#   - bash 4.0+ (for arrays and ${!var} syntax)
#   - Standard POSIX utilities (sed, xargs, mktemp)
#
# Author: Vatsal Sanjay
# Organization: CoMPhy Lab, Durham University

# Global arrays populated by extract_sweep_variables
SWEEP_VARS=()
SWEEP_VALUES=()

# Global variables populated by setup_sweep_temp_dir and generate_sweep_combinations
SWEEP_TEMP_DIR=""
SWEEP_CASE_NUM=0
SWEEP_COMBINATION_COUNT=0
SWEEP_CASE_FILES=()

# ============================================================
# Extract SWEEP_* variables from a sweep file
# Sets: SWEEP_VARS[], SWEEP_VALUES[]
# Usage: extract_sweep_variables <sweep_file>
# ============================================================
extract_sweep_variables() {
    local sweep_file="$1"

    if [ ! -f "$sweep_file" ]; then
        echo "ERROR: Sweep file not found: $sweep_file" >&2
        return 1
    fi

    # Reset arrays
    SWEEP_VARS=()
    SWEEP_VALUES=()

    # Read sweep file and extract SWEEP_* variables
    while IFS='=' read -r key value; do
        # Skip comments and empty lines
        [[ "${key:-}" =~ ^[[:space:]]*# ]] && continue
        [[ -z "${key:-}" ]] && continue

        # Match SWEEP_* variables
        if [[ "$key" =~ ^[[:space:]]*SWEEP_([^=]+) ]]; then
            local var_name="${BASH_REMATCH[1]}"
            # Remove inline comments and whitespace
            value=$(echo "${value:-}" | sed 's/#.*//' | xargs)

            SWEEP_VARS+=("$var_name")
            SWEEP_VALUES+=("$value")
        fi
    done < "$sweep_file"

    if [ ${#SWEEP_VARS[@]} -eq 0 ]; then
        echo "ERROR: No SWEEP_* variables found in $sweep_file" >&2
        return 1
    fi

    return 0
}

# ============================================================
# Validate case number range (4-digit, start <= end)
# Usage: validate_case_range <start> <end>
# ============================================================
validate_case_range() {
    local case_start="$1"
    local case_end="$2"

    if [ -z "$case_start" ] || [ -z "$case_end" ]; then
        echo "ERROR: CASE_START and CASE_END must be defined" >&2
        return 1
    fi

    if [ "$case_start" -lt 1000 ] || [ "$case_start" -gt 9999 ]; then
        echo "ERROR: CASE_START must be 4-digit (1000-9999), got: $case_start" >&2
        return 1
    fi

    if [ "$case_end" -lt "$case_start" ] || [ "$case_end" -gt 9999 ]; then
        echo "ERROR: CASE_END must be >= CASE_START and <= 9999, got: $case_end" >&2
        return 1
    fi

    return 0
}

# ============================================================
# Print sweep variables (for debugging/logging)
# Usage: print_sweep_variables
# ============================================================
print_sweep_variables() {
    echo "Sweep variables:"
    for i in "${!SWEEP_VARS[@]}"; do
        echo "  ${SWEEP_VARS[$i]} = ${SWEEP_VALUES[$i]}"
    done
}

# ============================================================
# Setup temporary directory for sweep generation
# Sets: SWEEP_TEMP_DIR
# Usage: setup_sweep_temp_dir [base_dir]
# ============================================================
setup_sweep_temp_dir() {
    local base_dir="${1:-}"

    if [ -n "$base_dir" ]; then
        # Use work directory (for HPC - avoids slow /tmp)
        SWEEP_TEMP_DIR="${base_dir}/.sweep_tmp_$$"
        mkdir -p "$SWEEP_TEMP_DIR" || {
            echo "ERROR: Failed to create temp directory: $SWEEP_TEMP_DIR" >&2
            return 1
        }
    else
        # Use system temp directory
        SWEEP_TEMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/sweep.XXXXXX")
    fi

    # Setup cleanup trap
    trap 'rm -rf "$SWEEP_TEMP_DIR"' EXIT

    return 0
}

# ============================================================
# Generate all parameter combinations (Cartesian product)
# Requires: SWEEP_VARS[], SWEEP_VALUES[], BASE_CONFIG, CASE_START
# Sets: SWEEP_CASE_FILES[], SWEEP_COMBINATION_COUNT
# Usage: generate_sweep_combinations [verbose]
# ============================================================
generate_sweep_combinations() {
    local verbose="${1:-0}"

    if [ -z "${SWEEP_TEMP_DIR:-}" ]; then
        echo "ERROR: Call setup_sweep_temp_dir first" >&2
        return 1
    fi

    if [ -z "${BASE_CONFIG:-}" ]; then
        echo "ERROR: BASE_CONFIG not defined" >&2
        return 1
    fi

    if [ ! -f "$BASE_CONFIG" ]; then
        echo "ERROR: Base configuration file not found: $BASE_CONFIG" >&2
        return 1
    fi

    # Initialize counters
    SWEEP_CASE_NUM=${CASE_START:-1000}
    SWEEP_COMBINATION_COUNT=0
    SWEEP_CASE_FILES=()

    # Recursive function to generate all combinations
    _generate_combinations_recursive() {
        local depth=$1
        shift
        local current_values=("$@")

        if [ $depth -eq ${#SWEEP_VARS[@]} ]; then
            # Base case: all variables assigned, create parameter file
            local case_file="${SWEEP_TEMP_DIR}/case_$(printf "%04d" "$SWEEP_CASE_NUM").params"

            # Copy base config
            cp "$BASE_CONFIG" "$case_file"

            # Override CaseNo
            if grep -q "^CaseNo=" "$case_file"; then
                sed -i'.bak' "s|^CaseNo=.*|CaseNo=${SWEEP_CASE_NUM}|" "$case_file"
            else
                echo "CaseNo=${SWEEP_CASE_NUM}" >> "$case_file"
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

            SWEEP_CASE_FILES+=("$case_file")
            SWEEP_COMBINATION_COUNT=$((SWEEP_COMBINATION_COUNT + 1))

            # Print summary if verbose
            if [ "$verbose" -eq 1 ]; then
                echo "Case $SWEEP_CASE_NUM:"
                for i in "${!SWEEP_VARS[@]}"; do
                    echo "  ${SWEEP_VARS[$i]} = ${current_values[$i]}"
                done
                echo ""
            fi

            SWEEP_CASE_NUM=$((SWEEP_CASE_NUM + 1))
            return
        fi

        # Recursive case: iterate through values for current variable
        local values="${SWEEP_VALUES[$depth]}"
        IFS=',' read -r -a value_array <<< "$values"

        for val in "${value_array[@]}"; do
            val=$(echo "$val" | xargs)  # Trim whitespace
            _generate_combinations_recursive $((depth + 1)) "${current_values[@]}" "$val"
        done
    }

    # Start recursion
    _generate_combinations_recursive 0

    return 0
}

# ============================================================
# Validate combination count matches expected range
# Usage: validate_combination_count <case_start> <case_end>
# ============================================================
validate_combination_count() {
    local case_start="$1"
    local case_end="$2"
    local expected_count=$((case_end - case_start + 1))

    if [ "$SWEEP_COMBINATION_COUNT" -ne "$expected_count" ]; then
        echo "WARNING: Generated $SWEEP_COMBINATION_COUNT combinations, but CASE_END suggests $expected_count" >&2
        echo "         Consider adjusting CASE_END in sweep file" >&2
    fi

    if [ "$SWEEP_COMBINATION_COUNT" -gt "$expected_count" ]; then
        echo "ERROR: Too many combinations ($SWEEP_COMBINATION_COUNT) for range $case_start-$case_end" >&2
        return 1
    fi

    return 0
}

# ============================================================
# Progress Tracking for Sweeps
# ============================================================

# Global variables for progress tracking
SWEEP_PROGRESS_START_TIME=0
SWEEP_PROGRESS_COMPLETED=0
SWEEP_PROGRESS_TOTAL=0

# Initialize progress tracking
# Usage: sweep_progress_init <total_cases>
sweep_progress_init() {
    SWEEP_PROGRESS_TOTAL="$1"
    SWEEP_PROGRESS_COMPLETED=0
    SWEEP_PROGRESS_START_TIME=$(date +%s)
}

# Update and display progress
# Usage: sweep_progress_update <case_no>
# Returns: Prints progress line with ETA
sweep_progress_update() {
    local case_no="$1"
    SWEEP_PROGRESS_COMPLETED=$((SWEEP_PROGRESS_COMPLETED + 1))

    local elapsed=$(($(date +%s) - SWEEP_PROGRESS_START_TIME))

    if [ $SWEEP_PROGRESS_COMPLETED -gt 1 ] && [ $elapsed -gt 0 ]; then
        local avg_time=$((elapsed / (SWEEP_PROGRESS_COMPLETED - 1)))
        local remaining=$((SWEEP_PROGRESS_TOTAL - SWEEP_PROGRESS_COMPLETED))
        local eta=$((avg_time * remaining))

        # Format ETA as HH:MM:SS
        local eta_h=$((eta / 3600))
        local eta_m=$(((eta % 3600) / 60))
        local eta_s=$((eta % 60))

        printf "Progress: %d/%d cases completed (ETA: %02d:%02d:%02d)\n" \
            "$SWEEP_PROGRESS_COMPLETED" "$SWEEP_PROGRESS_TOTAL" \
            "$eta_h" "$eta_m" "$eta_s"
    else
        printf "Progress: %d/%d cases completed\n" \
            "$SWEEP_PROGRESS_COMPLETED" "$SWEEP_PROGRESS_TOTAL"
    fi
}
