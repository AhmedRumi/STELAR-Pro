#!/usr/bin/env bash
# run-bulk-simulated.sh
#
# Runs sim.sh and test-stelarx-simulated.sh or test-baseline-simulated.sh
# over all combinations of parameter lists.
#
# Usage:
#   ./run-bulk-simulated.sh -m stelar
#   ./run-bulk-simulated.sh --project-root /path/to/checkout

set -euo pipefail

STELARX_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${STELARX_ROOT}/scripts/phylogeny-data-dir.sh"

BASE_DIR=""
BASE_DIR_PROVIDED=false
METHOD="stelarx"  # default method
FRESH=false
NO_NOTIFY=false
GPU_MONITOR=true
SIMPHY_DATA_DIR=""
NUM_REPLICATES=1

T_LIST=(10)

# T_LIST=(1000 2500 5000 7500 10000)

G_LIST=(1000)
SB_LIST=(0.000001)
SPMIN_LIST=(100000)
SPMAX_LIST=(200000)

# Method-specific options (passed through)
ASTER_OPTS=""
ASTER_BIN=""
ASTRAL_OPTS=""
STELARX_OPTS_LIST_RAW=""
ASTRAL_XMS=""
ASTRAL_XMX=""
TREEQMC_OPTS=""
WQFM_OPTS=""
SUPERTRIPLETS_OPTS=""
TMC_OPTS=""

print_help() {
  cat <<EOF
run-bulk-simulated.sh

Runs sim.sh and test-stelarx-simulated.sh or test-baseline-simulated.sh for all combinations of parameter lists.

Options:
  --method, -m      Method to use: stelarx (default: stelarx)
  --project-root    STELAR-X checkout root (default: this script's directory)
  --base-dir, -b    Compatibility alias for --project-root
  --num-replicates, -n  Number of replicates to run (default: 1)
  --taxa-list LIST       Comma/space-separated taxon counts (default: 10)
  --genes-list LIST      Comma/space-separated gene-tree counts (default: 10)
  --sb-list LIST         Comma/space-separated speciation rates
  --spmin-list LIST      Comma/space-separated minimum population sizes
  --spmax-list LIST      Comma/space-separated maximum population sizes
  --simphy-data-dir DIR  Store/read generated datasets under DIR
                         (default: \$PHYLOGENY_DATA_DIR/simphy/data)
  --fresh           Pass --fresh to sim.sh and test scripts (recreate outputs)
  --no-gpu-monitor  Disable GPU-memory sampling
  --no-notify, -nn  Disable completion notifications
  --opts, --alg-opts       Extra options for one STELAR-X simulated setting
  --opts-list, --alg-opts-list
                         Semicolon-separated list of STELAR-X option strings to loop over
  --help, -h        Show this message

Examples:
  ./run-bulk-simulated.sh --opts "--search-space S2 -vv"
  ./run-bulk-simulated.sh --taxa-list "10,20" --genes-list "10,50" --num-replicates 3
  ./run-bulk-simulated.sh --opts-list "--search-space S1 -vv;--search-space S2 -vv;--search-space S3 -vv"
EOF
}

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --method|-m) METHOD="$2"; shift 2 ;;
    --project-root|--base-dir|-b) BASE_DIR="$2"; BASE_DIR_PROVIDED=true; shift 2 ;;
    --num-replicates|-n) NUM_REPLICATES="$2"; shift 2 ;;
    --taxa-list) read -r -a T_LIST <<< "${2//,/ }"; shift 2 ;;
    --genes-list) read -r -a G_LIST <<< "${2//,/ }"; shift 2 ;;
    --sb-list) read -r -a SB_LIST <<< "${2//,/ }"; shift 2 ;;
    --spmin-list) read -r -a SPMIN_LIST <<< "${2//,/ }"; shift 2 ;;
    --spmax-list) read -r -a SPMAX_LIST <<< "${2//,/ }"; shift 2 ;;
    --simphy-data-dir) SIMPHY_DATA_DIR="$2"; shift 2 ;;
    --opts|--alg-opts|--stelarx-opts) ASTRAL_OPTS="$2"; shift 2 ;;
    --opts=*|--alg-opts=*|--stelarx-opts=*) ASTRAL_OPTS="${1#*=}"; shift ;;
    --opts-list|--alg-opts-list|--stelarx-opts-list) STELARX_OPTS_LIST_RAW="$2"; shift 2 ;;
    --opts-list=*|--alg-opts-list=*|--stelarx-opts-list=*) STELARX_OPTS_LIST_RAW="${1#*=}"; shift ;;
    --fresh) FRESH=true; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

# Validate method
case "$METHOD" in
  stelarx|astral-x|stelar) METHOD="stelarx" ;;
  *)
    echo "Error: --method must be stelarx."
    exit 1
    ;;
esac

# Keep the defaults deliberately small. Larger experiment matrices must be
# requested explicitly through the list options above.
for parameter_list in T_LIST G_LIST SB_LIST SPMIN_LIST SPMAX_LIST; do
  declare -n values="$parameter_list"
  if [[ ${#values[@]} -eq 0 ]]; then
    echo "Error: $parameter_list cannot be empty."
    exit 1
  fi
done
unset -n values

if [[ ! "$NUM_REPLICATES" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: --num-replicates must be a positive integer."
  exit 1
fi

SIMPHY_DATA_DIR="$(stelarx_prepare_simphy_data_dir "$SIMPHY_DATA_DIR")"

# -------------------------------
# execution
# -------------------------------

# Build base-dir argument if provided
BASE_DIR_ARGS=()
if $BASE_DIR_PROVIDED; then
  BASE_DIR_ARGS=(--project-root "$BASE_DIR")
  echo "Project root: $BASE_DIR"
else
  echo "Project root: $STELARX_ROOT"
fi

# Build fresh argument if provided
FRESH_ARGS=()
if $FRESH; then
  FRESH_ARGS=(--fresh)
  echo "Fresh:    yes"
else
  echo "Fresh:    no"
fi
SIM_DATA_ARGS=(--simphy-data-dir "$SIMPHY_DATA_DIR")
SHARED_TEST_ARGS=("${SIM_DATA_ARGS[@]}")
if [[ "$GPU_MONITOR" == false ]]; then
  SHARED_TEST_ARGS+=(--no-gpu-monitor)
fi
if [[ "$NO_NOTIFY" == true ]]; then
  SHARED_TEST_ARGS+=(--no-notify)
fi
echo "Method:   $METHOD"
echo "Replicates: $NUM_REPLICATES"
echo "SimPhy data: $SIMPHY_DATA_DIR"

STELARX_OPTS_LIST=()
if [[ -n "$STELARX_OPTS_LIST_RAW" ]]; then
  IFS=';' read -r -a raw_opts_list <<< "$STELARX_OPTS_LIST_RAW"
  for opts in "${raw_opts_list[@]}"; do
    opts="$(echo "$opts" | sed 's/^ *//;s/ *$//')"
    [[ -n "$opts" ]] && STELARX_OPTS_LIST+=("$opts")
  done
fi
if [[ ${#STELARX_OPTS_LIST[@]} -eq 0 ]]; then
  STELARX_OPTS_LIST+=("${ASTRAL_OPTS}")
fi

echo "Starting bulk runs..."

for t in "${T_LIST[@]}"; do
  for g in "${G_LIST[@]}"; do
    for sb in "${SB_LIST[@]}"; do
      for spmin in "${SPMIN_LIST[@]}"; do
        for spmax in "${SPMAX_LIST[@]}"; do

          echo ">>> Running: t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax (method=$METHOD)"
          
          ./sim.sh -rs "$NUM_REPLICATES" "${BASE_DIR_ARGS[@]}" "${SIM_DATA_ARGS[@]}" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${FRESH_ARGS[@]}"
          
          # Run replicates
          for ((i=1; i<=NUM_REPLICATES; i++)); do
            echo "  Running replicate R$i with $METHOD"
            
            for STELARX_OPTS_ITEM in "${STELARX_OPTS_LIST[@]}"; do
              TEST_CMD=("${STELARX_ROOT}/test-stelarx-simulated.sh" -r "R$i" "${BASE_DIR_ARGS[@]}" "${SHARED_TEST_ARGS[@]}" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${FRESH_ARGS[@]}")
              if [[ -n "$STELARX_OPTS_ITEM" ]]; then
                TEST_CMD+=(--opts "$STELARX_OPTS_ITEM")
              fi
              "${TEST_CMD[@]}"
            done
          done

        done
      done
    done
  done
done

echo "All runs finished."














T_LIST=(1000)
G_LIST=(1000 2500 5000 7500 10000)
SB_LIST=(0.000001)
SPMIN_LIST=(100000)
SPMAX_LIST=(200000)












echo "Starting bulk runs... phase 2"

for t in "${T_LIST[@]}"; do
  for g in "${G_LIST[@]}"; do
    for sb in "${SB_LIST[@]}"; do
      for spmin in "${SPMIN_LIST[@]}"; do
        for spmax in "${SPMAX_LIST[@]}"; do

          echo ">>> Running: t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax (method=$METHOD)"
          
          ./sim.sh -rs "$NUM_REPLICATES" "${BASE_DIR_ARGS[@]}" "${SIM_DATA_ARGS[@]}" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${FRESH_ARGS[@]}"
          
          # Run replicates
          for ((i=1; i<=NUM_REPLICATES; i++)); do
            echo "  Running replicate R$i with $METHOD"
            
            for STELARX_OPTS_ITEM in "${STELARX_OPTS_LIST[@]}"; do
              TEST_CMD=("${STELARX_ROOT}/test-stelarx-simulated.sh" -r "R$i" "${BASE_DIR_ARGS[@]}" "${SHARED_TEST_ARGS[@]}" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${FRESH_ARGS[@]}")
              if [[ -n "$STELARX_OPTS_ITEM" ]]; then
                TEST_CMD+=(--opts "$STELARX_OPTS_ITEM")
              fi
              "${TEST_CMD[@]}"
            done
          done

        done
      done
    done
  done
done

echo "All runs finished."
