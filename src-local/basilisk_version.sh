#!/bin/bash
# basilisk_version.sh - Centralized Basilisk version configuration
#
# Description:
#   Defines the Basilisk version and provides a helper function for installation.
#   Source this file in all scripts that need to install or reference Basilisk.
#
# Usage:
#   source src-local/basilisk_version.sh
#   install_basilisk              # Install with default settings
#   install_basilisk --hard       # Force reinstall
#
# Author: Vatsal Sanjay
# Organization: CoMPhy Lab, Durham University

# ============================================================
# Basilisk Version Configuration
# ============================================================
# Pin to a specific version for reproducibility
# Update this value when upgrading Basilisk across the project
BASILISK_REF="v2026-01-13"

# Installation script URL (from comphy-lab/basilisk-C)
BASILISK_INSTALL_URL="https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh"

# ============================================================
# Installation Helper
# ============================================================
# Install Basilisk with the pinned version
# Usage: install_basilisk [options]
# Options are passed directly to the install script (e.g., --hard)
install_basilisk() {
    if ! command -v curl >/dev/null 2>&1; then
        echo "ERROR: curl not found (required to install Basilisk)" >&2
        return 1
    fi

    echo "Installing Basilisk (ref: ${BASILISK_REF})..."
    curl -sL "$BASILISK_INSTALL_URL" | bash -s -- --ref="$BASILISK_REF" "$@"
}
