#!/usr/bin/env bash
# Compatibility entry point for the current rooted-triplet validation suite.
# The pre-migration runner used an unrooted quartet oracle and accepted fixtures
# that STELAR-X now intentionally rejects, so keeping that logic here produced
# false failures against the current command contract.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
ARGS=()

usage() {
  cat <<'EOF'
Usage: test/run_tests.sh [--cpu|--gpu] [--quick] [--skip-packaging]

Compatibility wrapper around run_stelarx_comprehensive_tests.sh.
  --cpu               Run CPU validation only.
  --gpu               Require the strict CUDA validation layer.
  --quick             Use the comprehensive suite's reduced randomized matrix.
  --skip-packaging    Skip the portable-application build and smoke test.

The old fixture filter and --search-mode options belonged to the superseded
unrooted-quartet runner and are no longer supported.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cpu) ARGS+=(--cpu-only); shift ;;
    --gpu) ARGS+=(--require-gpu); shift ;;
    --quick|--skip-packaging|--cpu-only|--require-gpu)
      ARGS+=("$1"); shift ;;
    --no-build)
      echo "Note: --no-build is ignored; comprehensive validation always rebuilds Java." >&2
      shift ;;
    -h|--help) usage; exit 0 ;;
    --search-mode)
      echo "Error: --search-mode belonged to the retired unrooted-quartet runner." >&2
      usage >&2
      exit 2 ;;
    --*)
      echo "Error: unknown option: $1" >&2
      usage >&2
      exit 2 ;;
    *)
      echo "Error: fixture filters belonged to the retired unrooted-quartet runner: $1" >&2
      usage >&2
      exit 2 ;;
  esac
done

exec "${ROOT}/test/run_stelarx_comprehensive_tests.sh" "${ARGS[@]}"
