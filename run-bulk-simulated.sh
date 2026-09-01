#!/usr/bin/env bash
# run-bulk-simulated.sh
#
# Runs sim.sh and test-stelar-pro-simulated.sh or test-baseline-simulated.sh
# over all combinations of parameter lists.
#
# Usage:
#   ./run-bulk-simulated.sh -m stelar
#   ./run-bulk-simulated.sh --project-root /path/to/checkout

set -euo pipefail

STELAR_PRO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${STELAR_PRO_ROOT}/scripts/phylogeny-data-dir.sh"

BASE_DIR=""
BASE_DIR_PROVIDED=false
METHOD="stelar-pro"  # default method
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

# Configurations that must not run through STELAR-Pro. Each entry is an exact
# TAXA,GENE_TREES,SB,SPMIN,SPMAX,REPLICATE tuple. Keep replicate names in R<n>
# form. sim.sh may still prepare the surrounding dataset batch; only the listed
# per-replicate inference/result run is purposefully skipped.
EXCLUDED_SIMULATED_CONFIGS=(
  "30000,1000,0.000001,100000,150000,R5"
  "30000,1000,0.000001,100000,250000,R4"
  "30000,1000,0.000001,100000,300000,R3"
  "40000,1000,0.000001,100000,150000,R2"
  "40000,1000,0.000001,100000,150000,R5"
  "40000,1000,0.000001,100000,200000,R5"
  "40000,1000,0.000001,100000,300000,R5"
  "1000,25000,0.000001,100000,200000,R3"
  "1000,25000,0.000001,100000,250000,R2"
  "75000,1000,0.000001,100000,200000,R5"
  "100000,1000,0.000001,100000,200000,R3"
  "100000,1000,0.000001,100000,200000,R4"
  "125000,1000,0.000001,100000,200000,R2"
  "125000,1000,0.000001,100000,200000,R3"
  "125000,1000,0.000001,100000,200000,R4"
  "150000,1000,0.000001,100000,200000,R4"
  "175000,1000,0.000001,100000,200000,R2"
  "175000,1000,0.000001,100000,200000,R4"
  "200000,1000,0.000001,100000,200000,R1"
  "200000,1000,0.000001,100000,200000,R2"
  "225000,1000,0.000001,100000,200000,R1"
  "225000,1000,0.000001,100000,200000,R3"
  "1000,75000,0.000001,100000,200000,R3"
  "1000,100000,0.000001,100000,200000,R4"
  "1000,125000,0.000001,100000,200000,R4"
  "1000,150000,0.000001,100000,200000,R2"
  "1000,200000,0.000001,100000,200000,R4"
  "1000,225000,0.000001,100000,200000,R2"
  "1000,275000,0.000001,100000,200000,R3"
  "1000,300000,0.000001,100000,200000,R2"
)

IS_SIMULATED_CONFIG_EXCLUDED() {
  local CANDIDATE_CONFIG="$1,$2,$3,$4,$5,$6"
  local EXCLUDED_CONFIG
  for EXCLUDED_CONFIG in "${EXCLUDED_SIMULATED_CONFIGS[@]}"; do
    [[ "$CANDIDATE_CONFIG" == "$EXCLUDED_CONFIG" ]] && return 0
  done
  return 1
}

# Method-specific options (passed through)
ASTER_OPTS=""
ASTER_BIN=""
ASTRAL_OPTS=""
STELAR_PRO_OPTS_LIST_RAW=""
ASTRAL_XMS=""
ASTRAL_XMX=""
TREEQMC_OPTS=""
WQFM_OPTS=""
SUPERTRIPLETS_OPTS=""
TMC_OPTS=""

# Permit the exclusion predicate and uppercase configuration array to be loaded
# by the isolated regression test without executing a simulated-data sweep.
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  return 0
fi

print_help() {
  cat <<EOF
run-bulk-simulated.sh

Runs sim.sh and test-stelar-pro-simulated.sh or test-baseline-simulated.sh for all combinations of parameter lists.

Options:
  --method, -m      Method to use: stelar-pro (default: stelar-pro)
  --project-root    STELAR-Pro checkout root (default: this script's directory)
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
  --opts, --alg-opts       Extra options for one STELAR-Pro simulated setting
  --opts-list, --alg-opts-list
                         Semicolon-separated list of STELAR-Pro option strings to loop over
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
    --opts|--alg-opts|--stelar-pro-opts) ASTRAL_OPTS="$2"; shift 2 ;;
    --opts=*|--alg-opts=*|--stelar-pro-opts=*) ASTRAL_OPTS="${1#*=}"; shift ;;
    --opts-list|--alg-opts-list|--stelar-pro-opts-list) STELAR_PRO_OPTS_LIST_RAW="$2"; shift 2 ;;
    --opts-list=*|--alg-opts-list=*|--stelar-pro-opts-list=*) STELAR_PRO_OPTS_LIST_RAW="${1#*=}"; shift ;;
    --fresh) FRESH=true; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

# Validate method
case "$METHOD" in
  stelar-pro|astral-x|stelar) METHOD="stelar-pro" ;;
  *)
    echo "Error: --method must be stelar-pro."
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

SIMPHY_DATA_DIR="$(stelar_pro_prepare_simphy_data_dir "$SIMPHY_DATA_DIR")"

# -------------------------------
# execution
# -------------------------------

# Build base-dir argument if provided
BASE_DIR_ARGS=()
if $BASE_DIR_PROVIDED; then
  BASE_DIR_ARGS=(--project-root "$BASE_DIR")
  echo "Project root: $BASE_DIR"
else
  echo "Project root: $STELAR_PRO_ROOT"
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

STELAR_PRO_OPTS_LIST=()
if [[ -n "$STELAR_PRO_OPTS_LIST_RAW" ]]; then
  IFS=';' read -r -a raw_opts_list <<< "$STELAR_PRO_OPTS_LIST_RAW"
  for opts in "${raw_opts_list[@]}"; do
    opts="$(echo "$opts" | sed 's/^ *//;s/ *$//')"
    [[ -n "$opts" ]] && STELAR_PRO_OPTS_LIST+=("$opts")
  done
fi
if [[ ${#STELAR_PRO_OPTS_LIST[@]} -eq 0 ]]; then
  STELAR_PRO_OPTS_LIST+=("${ASTRAL_OPTS}")
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
            REPLICATE_NAME="R$i"
            if IS_SIMULATED_CONFIG_EXCLUDED \
                "$t" "$g" "$sb" "$spmin" "$spmax" "$REPLICATE_NAME"; then
              echo "  SKIPPING excluded configuration: t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax replicate=$REPLICATE_NAME"
              continue
            fi
            echo "  Running replicate $REPLICATE_NAME with $METHOD"
            
            for STELAR_PRO_OPTS_ITEM in "${STELAR_PRO_OPTS_LIST[@]}"; do
              TEST_CMD=("${STELAR_PRO_ROOT}/test-stelar-pro-simulated.sh" -r "$REPLICATE_NAME" "${BASE_DIR_ARGS[@]}" "${SHARED_TEST_ARGS[@]}" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${FRESH_ARGS[@]}")
              if [[ -n "$STELAR_PRO_OPTS_ITEM" ]]; then
                TEST_CMD+=(--opts "$STELAR_PRO_OPTS_ITEM")
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
            REPLICATE_NAME="R$i"
            if IS_SIMULATED_CONFIG_EXCLUDED \
                "$t" "$g" "$sb" "$spmin" "$spmax" "$REPLICATE_NAME"; then
              echo "  SKIPPING excluded configuration: t=$t g=$g sb=$sb spmin=$spmin spmax=$spmax replicate=$REPLICATE_NAME"
              continue
            fi
            echo "  Running replicate $REPLICATE_NAME with $METHOD"
            
            for STELAR_PRO_OPTS_ITEM in "${STELAR_PRO_OPTS_LIST[@]}"; do
              TEST_CMD=("${STELAR_PRO_ROOT}/test-stelar-pro-simulated.sh" -r "$REPLICATE_NAME" "${BASE_DIR_ARGS[@]}" "${SHARED_TEST_ARGS[@]}" -t "$t" -g "$g" --sb "$sb" --spmin "$spmin" --spmax "$spmax" "${FRESH_ARGS[@]}")
              if [[ -n "$STELAR_PRO_OPTS_ITEM" ]]; then
                TEST_CMD+=(--opts "$STELAR_PRO_OPTS_ITEM")
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
