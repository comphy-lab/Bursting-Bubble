#!/bin/bash
# parse_params.sh - Shell library for parameter file parsing
#
# Description:
#   Provides reusable functions for parsing key=value parameter files,
#   validating required parameters, and generating parameter sweep cases.
#   This library is sourced by simulation scripts to handle configuration.
#
# Functions:
#   parse_param_file <file>     - Parse file, export as PARAM_* env vars
#   get_param <key> [default]   - Get parameter value with optional default
#   validate_required_params    - Check required parameters are set
#   print_params                - Debug: print all loaded parameters
#
# Usage:
#   source src-local/parse_params.sh
#   parse_param_file "config.params"
#   Oh=$(get_param "Oh" "0.01")
#
# Dependencies:
#   - bash 4.0+ (for associative arrays and ${!var} syntax)
#   - Standard POSIX utilities (sed, xargs, mktemp)
#
# Author: Vatsal Sanjay
# Organization: CoMPhy Lab, Durham University

# Parse a parameter file and export all parameters as environment variables
# Usage: parse_param_file <file>
parse_param_file() {
    local param_file=$1

    if [ ! -f "$param_file" ]; then
        echo "ERROR: Parameter file $param_file not found" >&2
        return 1
    fi

    # Read parameters (skip comments and empty lines)
    while IFS='=' read -r key value || [ -n "$key" ]; do
        # Skip comments and empty lines
        [[ "$key" =~ ^[[:space:]]*# ]] && continue
        [[ -z "$key" ]] && continue

        # Remove inline comments and whitespace
        value=$(echo "$value" | sed 's/#.*//' | xargs)
        key=$(echo "$key" | xargs)

        # Skip if key or value is empty
        [ -z "$key" ] && continue
        [ -z "$value" ] && continue

        # Export as environment variable with PARAM_ prefix
        export "PARAM_${key}=${value}"
    done < "$param_file"

    return 0
}

# Get a parameter value with optional default
# Usage: get_param <key> [default]
get_param() {
    local key=$1
    local default=${2:-}
    local var_name="PARAM_${key}"
    echo "${!var_name:-$default}"
}

# Validate that required variables are set in parameter file
# Usage: validate_required_params <param1> <param2> ...
validate_required_params() {
    local missing=0

    for param in "$@"; do
        local var_name="PARAM_${param}"
        if [ -z "${!var_name}" ]; then
            echo "ERROR: Required parameter '$param' not found" >&2
            missing=1
        fi
    done

    return $missing
}

# Print all loaded parameters (for debugging)
print_params() {
    echo "Loaded parameters:"
    env | grep "^PARAM_" | sort | while IFS='=' read -r key value; do
        key="${key#PARAM_}"
        echo "  $key = $value"
    done
}

# Validate a restart file exists, is non-empty, and is readable
# Usage: validate_restart_file [file_path]
# Returns 0 if valid, 1 if invalid
validate_restart_file() {
    local restart_file="${1:-restart}"

    if [ ! -f "$restart_file" ]; then
        echo "ERROR: Restart file not found: $restart_file" >&2
        return 1
    fi

    if [ ! -s "$restart_file" ]; then
        echo "ERROR: Restart file is empty: $restart_file" >&2
        return 1
    fi

    if [ ! -r "$restart_file" ]; then
        echo "ERROR: Restart file not readable: $restart_file" >&2
        return 1
    fi

    return 0
}

# Resolve Solver= from the parsed params into SOLVER_STEM / SOURCE_FILE_NAME /
# EXECUTABLE_NAME. Known stems must match files in simulationCases/.
# Usage: resolve_solver
resolve_solver() {
    local solver
    solver=$(get_param "Solver" "burstingBubble")
    solver="${solver%.c}"
    case "$solver" in
        burstingBubble|burstingBubble-drillResolution|burstingBubbleVE|burstingBubbleVE-drillResolution)
            ;;
        *)
            echo "ERROR: Unknown Solver='$solver'" >&2
            echo "Allowed: burstingBubble, burstingBubble-drillResolution, burstingBubbleVE, burstingBubbleVE-drillResolution" >&2
            return 1
            ;;
    esac
    SOLVER_STEM="$solver"
    SOURCE_FILE_NAME="${solver}.c"
    EXECUTABLE_NAME="$solver"
}

# Drill solvers must not enable Basilisk's Linux FPE trap (-D_GNU_SOURCE).
# -D_DEFAULT_SOURCE still exposes madvise() without re-arming the trap.
solver_uses_drill() {
    case "${SOLVER_STEM:-}" in
        *drillResolution) return 0 ;;
        *) return 1 ;;
    esac
}

solver_mpi_cc99() {
    if solver_uses_drill; then
        echo "mpicc -std=c99 -D_DEFAULT_SOURCE"
    else
        echo "mpicc -std=c99 -D_GNU_SOURCE=1"
    fi
}

# Append host- and case-dependent -D flags for qcc. Call after the
# case.params file has been parsed. two-phase*.h use #ifdef FILTERED,
# so FILTERED=0 must leave the macro undefined.
# Usage: append_solver_qcc_flags
append_solver_qcc_flags() {
    QCC_FLAGS="${QCC_FLAGS:-}"
    if [ -n "${BASILISK:-}" ] && grep -q "void set_prolongation" \
        "${BASILISK}/grid/multigrid-common.h" 2>/dev/null; then
        QCC_FLAGS="${QCC_FLAGS} -DVE_USE_SET_PROLONGATION"
    fi
    local filtered
    filtered=$(get_param "FILTERED" "1")
    if [ "$filtered" != "0" ] && [ "$filtered" != "1" ]; then
        echo "ERROR: FILTERED must be 0 or 1, got: $filtered" >&2
        return 1
    fi
    # Incoming QCC_FLAGS may already contain -DFILTERED; FILTERED=0 must
    # drop it so #ifdef FILTERED stays false.
    QCC_FLAGS="${QCC_FLAGS//-DFILTERED/}"
    if [ "$filtered" = "1" ]; then
        QCC_FLAGS="${QCC_FLAGS} -DFILTERED"
    fi
    export QCC_FLAGS
}
