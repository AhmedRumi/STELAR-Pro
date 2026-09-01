#!/usr/bin/env bash
# run_comparison_tests.sh — Head-to-head comparison of STELAR-Pro vs ASTRAL-MP.
#
# For each test case it:
#   1. Runs STELAR-Pro   with --dump-clusters  → cluster dump + species tree
#   2. Runs ASTRAL-MP  with --dump-clusters -p 0  → cluster dump + species tree
#      (-p 0 disables UPGMA/greedy-consensus enrichment so X is gene-tree-only,
#       the same search space STELAR-Pro builds)
#   3. Calls compare_with_astralmp.py to diff the two cluster dumps and trees
#
# Usage:
#   bash test/run_comparison_tests.sh [TC_FILTER] [options]
#
# Options:
#   --cpu               force CPU mode for STELAR-Pro (default: --cpu for safety)
#   --gpu               use GPU for STELAR-Pro
#   --no-build          skip building both tools before running
#   --complete          add --autocomplete-incomplete-gene-trees to STELAR-Pro
#                       (use for TCs with missing taxa)
#
# TC_FILTER: glob pattern matching input file basenames, e.g. "tc1*" or "tc[13]"
#            default: all tc*.tre files in test/input/
#
# Exit code: 0 if all TCs pass, 1 if any fail.
#
# Requires an ASTRAL-MP checkout. Set ASTRALMP_ROOT to override the default
# <STELAR-Pro checkout>/astral-my location.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
INPUT_DIR="${SCRIPT_DIR}/input"
TMP_DIR="${SCRIPT_DIR}/output/comparison"

BUILD_DIR="${ROOT_DIR}/build"
NATIVE_DIR="${ROOT_DIR}/native"

ASTRALMP_ROOT="${ASTRALMP_ROOT:-${ROOT_DIR}/astral-my}"
ASTRALMP_RUN="${ASTRALMP_ROOT}/ASTRAL/run_astral.sh"

COMPARE_PY="${SCRIPT_DIR}/compare_with_astralmp.py"

# ── defaults ──────────────────────────────────────────────────────────────────
FILTER="tc[0-9]*"
COMPUTE_MODE="--cpu"
SKIP_BUILD=0
AUTOCOMPLETE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpu)       COMPUTE_MODE="--cpu"; shift ;;
        --gpu)       COMPUTE_MODE="--gpu"; shift ;;
        --no-build)  SKIP_BUILD=1; shift ;;
        --complete)  AUTOCOMPLETE="--autocomplete-incomplete-gene-trees"; shift ;;
        *)           FILTER="$1"; shift ;;
    esac
done

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m'

pass=0
fail=0
skip=0

mkdir -p "$TMP_DIR"

# ── build phase ───────────────────────────────────────────────────────────────
if [[ $SKIP_BUILD -eq 0 ]]; then
    echo -e "${BOLD}Building STELAR-Pro...${NC}"
    bash "${ROOT_DIR}/build.sh"
    bash "${ROOT_DIR}/build_native.sh" 2>/dev/null || true

    echo -e "${BOLD}Building ASTRAL-MP...${NC}"
    bash "${ASTRALMP_ROOT}/ASTRAL/compile_astral.sh"
    echo ""
fi

STELAR_PRO_CMD=(java -Djava.library.path="$NATIVE_DIR" -cp "$BUILD_DIR" stelarx.Main
             $COMPUTE_MODE --search-mode local)

# ── per-TC runner ─────────────────────────────────────────────────────────────
run_tc() {
    local input="$1"
    local name
    name="$(basename "$input" .tre)"

    # skip _true ground-truth files and non-TC files
    [[ "$name" == *_true ]] && return 0
    [[ "$name" == test_* ]] && return 0

    printf "  %-45s" "$name"

    local cx_file="${TMP_DIR}/${name}_stelarx_clusters.txt"
    local cm_file="${TMP_DIR}/${name}_astralmp_clusters.txt"
    local tx_file="${TMP_DIR}/${name}_stelarx_tree.tre"
    local tm_file="${TMP_DIR}/${name}_astralmp_tree.tre"

    # ── Run STELAR-Pro ──────────────────────────────────────────────────────────
    local stelarx_log
    if ! stelarx_log="$("${STELAR_PRO_CMD[@]}" \
            -i "$input" \
            -o "$tx_file" \
            --dump-clusters "$cx_file" \
            ${AUTOCOMPLETE} 2>&1)"; then
        printf "${RED}FAIL${NC}  STELAR-Pro error\n"
        echo "$stelarx_log" | tail -5 | sed 's/^/    /'
        ((fail++)) || true
        return 0
    fi

    # ── Run ASTRAL-MP (p=0: no extra bipartitions) ────────────────────────────
    local astralmp_log
    if ! astralmp_log="$(bash "$ASTRALMP_RUN" \
            -i "$(realpath "$input")" \
            -o "$(realpath "$tm_file")" \
            --dump-clusters "$(realpath "$cm_file")" \
            -p 0 -C -t 0 2>&1)"; then
        printf "${RED}FAIL${NC}  ASTRAL-MP error\n"
        echo "$astralmp_log" | tail -5 | sed 's/^/    /'
        ((fail++)) || true
        return 0
    fi

    # ── Compare ───────────────────────────────────────────────────────────────
    local cmp_out
    if cmp_out="$(python3 "$COMPARE_PY" \
            --clusters-stelar-pro  "$cx_file" \
            --clusters-astralmp "$cm_file" \
            --tree-stelar-pro      "$tx_file" \
            --tree-astralmp     "$tm_file" \
            --label             "$name" 2>&1)"; then
        printf "${GREEN}PASS${NC}\n"
        ((pass++)) || true
    else
        printf "${RED}FAIL${NC}\n"
        echo "$cmp_out" | sed 's/^/    /'
        ((fail++)) || true
    fi
}

# ── iterate inputs ────────────────────────────────────────────────────────────
echo -e "${BOLD}Running STELAR-Pro vs ASTRAL-MP comparison tests${NC}"
echo "  Input dir : $INPUT_DIR"
echo "  Tmp dir   : $TMP_DIR"
echo "  Mode      : STELAR-Pro $COMPUTE_MODE / ASTRAL-MP CPU-only -p 0"
[[ -n "$AUTOCOMPLETE" ]] && echo "  Autocomplete: on"
echo ""

shopt -s nullglob
# NOTE: FILTER must not be quoted in the glob so bash expands it
inputs=("${INPUT_DIR}"/${FILTER}.tre)
if [[ ${#inputs[@]} -eq 0 ]]; then
    echo "No input files matching '${FILTER}.tre' in $INPUT_DIR"
    exit 1
fi

for f in "${inputs[@]}"; do
    run_tc "$f"
done

echo ""
echo -e "${BOLD}Results: ${GREEN}${pass} passed${NC}  ${RED}${fail} failed${NC}"
[[ $fail -eq 0 ]] && exit 0 || exit 1
