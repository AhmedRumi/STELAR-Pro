#!/bin/bash
#
# STELAR-Pro Build Script
# ======================
# Quick rebuild of STELAR-Pro (Java + CUDA).
# Skips Java/JDK checks for fast iteration during development.
#
# Usage:
#   ./build.sh              # Quick rebuild (auto-detects CUDA)
#   ./build.sh --clean      # Clean rebuild from scratch
#   ./build.sh --no-cuda    # Skip CUDA compilation
#
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/install.sh" --quick "$@"
