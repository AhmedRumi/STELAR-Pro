#!/usr/bin/env bash
# Run ASTRAL-Pro3 on one repository SimPhy replicate and write the same
# per-replicate experiment CSV/checkpoint layout used by STELAR-Pro.

set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_ROOT}/scripts/phylogeny-data-dir.sh"
source "${SCRIPT_ROOT}/scripts/simulated-outputs-mirror.sh"
source "${SCRIPT_ROOT}/experiment-setting-name.sh"
ORIGINAL_INVOCATION=("$0" "$@")

NTFY_CHANNEL_NAME="${NTFY_CHANNEL_NAME:-anik-phylo-stx}"
TAXA_NUM=""
GENE_TREES=""
REPLICATE=R1
BASE_DIR="$SCRIPT_ROOT"
SIMPHY_DIR=""
SIMPHY_DIR_SET=false
SIMPHY_DATA_DIR=""
SIMULATED_OUTPUTS_DIR=""
NO_OUTPUTS_MIRROR=false
SB=0.000001
SPMIN=500000
SPMAX=1500000
USE_LEGACY_LAYOUT=false
INCOMPLETE=false
ASTRAL_PRO3_OPTS=""
ASTRAL_PRO3_BIN=""
FRESH=false
TIME_MONITOR=true
GPU_MONITOR=true
NO_NOTIFY=false
DEBUG=0

print_help() {
  cat <<EOF
test-astral-pro3-simulated.sh

Required:
  --taxa-num, -t N       Number of taxa
  --gene-trees, -g N     Number of gene trees

Optional:
  --replicate, -r NAME   Replicate (default: R1)
  --project-root DIR     STELAR-Pro checkout root
  --base-dir, -b DIR     Compatibility alias for --project-root
  --simphy-dir DIR       SimPhy installation directory
  --simphy-data-dir DIR  Simulated-data root
  --simulated-outputs-dir DIR
                         Results mirror root
                         (default: \$PHYLOGENY_DATA_DIR/outputs/gdl-simulation)
  --gdl-simulation-outputs-dir DIR
                         Alias for --simulated-outputs-dir
  --simphy-outputs-dir DIR
                         Compatibility alias for --simulated-outputs-dir
  --no-outputs-mirror   Do not copy this result into the outputs mirror
  --sb VALUE             Speciation/birthrate parameter
  --spmin VALUE          Minimum population size
  --spmax VALUE          Maximum population size
  --opts "..."           Extra ASTRAL-Pro3 options
  --astral-pro3-bin FILE Executable override
  --use-legacy-layout    Use TAXA_GENES/REPLICATE data layout
  --incomplete           Use the matching _incomplete dataset
  --fresh                Rerun a completed setting
  --no-time-monitor      Disable time/RAM monitoring
  --no-gpu-monitor       Disable GPU monitoring
  --no-notify, -nn       Disable completion notifications
  --debug                Enable shell tracing
  --help, -h             Show this message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --taxa-num|--taxa_num|-t) TAXA_NUM="$2"; shift 2 ;;
    --gene-trees|--gene_trees|-g) GENE_TREES="$2"; shift 2 ;;
    --replicate|-r) REPLICATE="$2"; shift 2 ;;
    --project-root|--base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --simphy-dir) SIMPHY_DIR="$2"; SIMPHY_DIR_SET=true; shift 2 ;;
    --simphy-data-dir) SIMPHY_DATA_DIR="$2"; shift 2 ;;
    --simulated-outputs-dir|--gdl-simulation-outputs-dir|--simphy-outputs-dir) SIMULATED_OUTPUTS_DIR="$2"; shift 2 ;;
    --no-outputs-mirror) NO_OUTPUTS_MIRROR=true; shift ;;
    --sb) SB="$2"; shift 2 ;;
    --spmin) SPMIN="$2"; shift 2 ;;
    --spmax) SPMAX="$2"; shift 2 ;;
    --opts|--alg-opts|--astral-pro3-opts) ASTRAL_PRO3_OPTS="$2"; shift 2 ;;
    --opts=*|--alg-opts=*|--astral-pro3-opts=*) ASTRAL_PRO3_OPTS="${1#*=}"; shift ;;
    --astral-pro3-bin|--astral-pro-bin) ASTRAL_PRO3_BIN="$2"; shift 2 ;;
    --use-legacy-layout) USE_LEGACY_LAYOUT=true; shift ;;
    --incomplete) INCOMPLETE=true; shift ;;
    --fresh) FRESH=true; shift ;;
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --debug) DEBUG=1; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Error: unknown option: $1" >&2; print_help >&2; exit 2 ;;
  esac
done

if [[ -z "$TAXA_NUM" || -z "$GENE_TREES" ]]; then
  echo "Error: --taxa-num and --gene-trees are required." >&2
  exit 2
fi
if [[ ! "$REPLICATE" =~ ^R[1-9][0-9]*$ ]]; then
  echo "Error: --replicate must have the form R1, R2, ..." >&2
  exit 2
fi

BASE_DIR="$(realpath "$BASE_DIR")"
[[ -n "$SIMPHY_DIR" ]] || SIMPHY_DIR="${BASE_DIR}/simphy"
SIMPHY_DIR="$(realpath -m "$SIMPHY_DIR")"
SIMPHY_DATA_DIR="$(stelar_pro_prepare_simphy_data_dir "$SIMPHY_DATA_DIR")"
if [[ "$NO_OUTPUTS_MIRROR" == false ]]; then
  SIMULATED_OUTPUTS_DIR="$(stelar_pro_simulated_outputs_dir "$SIMPHY_DATA_DIR" "$SIMULATED_OUTPUTS_DIR")"
fi
SETTING_NAME="$(build_astral_pro3_setting_name_from_opts "$ASTRAL_PRO3_OPTS")"

if [[ "$NO_OUTPUTS_MIRROR" == false ]]; then
  echo "Outputs mirror: $SIMULATED_OUTPUTS_DIR"
fi

if [[ "$USE_LEGACY_LAYOUT" == true ]]; then
  RUN_DIR="${SIMPHY_DATA_DIR}/${TAXA_NUM}_${GENE_TREES}/${REPLICATE}"
else
  DATASET_NAME="t_${TAXA_NUM}_g_${GENE_TREES}_sb_${SB}_spmin_${SPMIN}_spmax_${SPMAX}"
  [[ "$INCOMPLETE" == true ]] && DATASET_NAME+="_incomplete"
  RUN_DIR="${SIMPHY_DATA_DIR}/${DATASET_NAME}/${REPLICATE}"
fi

GENE_TREES_FILE="${RUN_DIR}/all_gt.tre"
TRUE_TREE="${RUN_DIR}/s_tree.trees"
RESULTS_DIR="${RUN_DIR}/astral-pro3-outputs/${SETTING_NAME}"
OUTPUT_TREE="${RESULTS_DIR}/out-astral-pro3.tre"
STAT_FILE="${RESULTS_DIR}/stat-astral-pro3.csv"
SUCCESS_FILE="${RESULTS_DIR}/.astral-pro3.success"
LOCK_FILE="${RESULTS_DIR}/.astral-pro3.lock"
WRAPPER_STATS="${OUTPUT_TREE%.tre}_stats.csv"
COMMAND_RECORD="${OUTPUT_TREE%.tre}.command"

mirror_results_leaf() {
  [[ "$NO_OUTPUTS_MIRROR" == true ]] && return 0
  if ! stelar_pro_mirror_simulated_results "$SIMPHY_DATA_DIR" "$SIMULATED_OUTPUTS_DIR" "$RESULTS_DIR"; then
    echo "WARNING: could not refresh outputs mirror for $RESULTS_DIR" >&2
  fi
  return 0
}

if [[ "$FRESH" == false && -s "$OUTPUT_TREE" && -f "$SUCCESS_FILE" && -f "$STAT_FILE" ]]; then
  mirror_results_leaf
  echo "SKIPPING: successful output already exists at $OUTPUT_TREE. Use --fresh to rerun."
  exit 0
fi

if [[ ! -f "$GENE_TREES_FILE" ]]; then
  if [[ "$USE_LEGACY_LAYOUT" == true ]]; then
    echo "Error: gene-tree file not found: $GENE_TREES_FILE" >&2
    exit 6
  fi

  REPLICATE_COUNT="${REPLICATE#R}"
  if [[ "$INCOMPLETE" == true ]]; then
    SIM_CMD=("${BASE_DIR}/sim_incomplete.sh" -t "$TAXA_NUM" -g "$GENE_TREES" -rs "$REPLICATE_COUNT"
      --sb "$SB" --spmin "$SPMIN" --spmax "$SPMAX" --simphy-data-dir "$SIMPHY_DATA_DIR")
    [[ "$SIMPHY_DIR_SET" == true ]] && SIM_CMD+=(--simphy-dir "$SIMPHY_DIR")
    [[ "$FRESH" == true ]] && SIM_CMD+=(--fresh-inc)
  else
    SIM_CMD=("${BASE_DIR}/sim.sh" -t "$TAXA_NUM" -g "$GENE_TREES" -r "$REPLICATE" -rs "$REPLICATE_COUNT"
      --sb "$SB" --spmin "$SPMIN" --spmax "$SPMAX" --simphy-data-dir "$SIMPHY_DATA_DIR")
    [[ "$SIMPHY_DIR_SET" == true ]] && SIM_CMD+=(--simphy-dir "$SIMPHY_DIR")
    [[ "$FRESH" == true ]] && SIM_CMD+=(--fresh)
  fi
  echo "==> Bootstrapping missing simulated dataset"
  "${SIM_CMD[@]}"
fi
if [[ ! -f "$GENE_TREES_FILE" ]]; then
  echo "Error: gene-tree file not found after dataset preparation: $GENE_TREES_FILE" >&2
  exit 6
fi
if [[ ! -f "$TRUE_TREE" ]]; then
  echo "Error: reference species tree not found: $TRUE_TREE" >&2
  exit 6
fi

mkdir -p "$RESULTS_DIR"
rm -f "$LOCK_FILE" "$SUCCESS_FILE"

CMD=("${BASE_DIR}/run-astral-pro3-with-monitor.sh"
  --input "$GENE_TREES_FILE" --output "$OUTPUT_TREE"
  --reference-species-tree "$TRUE_TREE" --no-notify)
[[ -n "$ASTRAL_PRO3_OPTS" ]] && CMD+=(--opts "$ASTRAL_PRO3_OPTS")
[[ -n "$ASTRAL_PRO3_BIN" ]] && CMD+=(--astral-pro3-bin "$ASTRAL_PRO3_BIN")
[[ "$TIME_MONITOR" == false ]] && CMD+=(--no-time-monitor)
[[ "$GPU_MONITOR" == false ]] && CMD+=(--no-gpu-monitor)
[[ "$DEBUG" == 1 ]] && CMD+=(--debug)

echo "==> Running ASTRAL-Pro3 on ${REPLICATE} [${SETTING_NAME}]"
set +e
"${CMD[@]}"
EXIT_CODE=$?
set -e

RUNNING_TIME=NA
MAX_CPU_MB=NA
MAX_GPU_MB=NA
RF_RATE=NA
if [[ -f "$WRAPPER_STATS" ]]; then
  RUNNING_TIME=$(awk -F, 'NR==2 {print $4}' "$WRAPPER_STATS")
  MAX_CPU_MB=$(awk -F, 'NR==2 {print $5}' "$WRAPPER_STATS")
  MAX_GPU_MB=$(awk -F, 'NR==2 {print $6}' "$WRAPPER_STATS")
  RF_RATE=$(awk -F, 'NR==2 {print $8}' "$WRAPPER_STATS")
fi

printf '%s\n' 'alg,setting,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,optimal-triplet-score,running-time-s,max-cpu-mb,max-gpu-mb' > "$STAT_FILE"
printf 'astral-pro3,%s,%s,%s,%s,%s,%s,%s,%s,NA,%s,%s,%s\n' \
  "$SETTING_NAME" "$TAXA_NUM" "$GENE_TREES" "$REPLICATE" "$SB" "$SPMIN" "$SPMAX" \
  "$RF_RATE" "$RUNNING_TIME" "$MAX_CPU_MB" "$MAX_GPU_MB" >> "$STAT_FILE"

if [[ -f "$COMMAND_RECORD" ]]; then
  {
    echo
    echo "# simulated run context"
    echo "# dataset: $(basename "$(dirname "$RUN_DIR")")"
    echo "# replicate: $REPLICATE"
    echo "# setting: $SETTING_NAME"
    echo "# true tree: $TRUE_TREE"
    echo "# RF rate: $RF_RATE"
    printf '# simulated wrapper invocation: '
    stelar_pro_print_shell_command "${ORIGINAL_INVOCATION[@]}"
    printf '# monitor wrapper command: '
    stelar_pro_print_shell_command "${CMD[@]}"
  } >> "$COMMAND_RECORD"
fi

if [[ $EXIT_CODE -eq 0 && -s "$OUTPUT_TREE" ]]; then
  touch "$LOCK_FILE"
  printf 'exit_code=0\noutput=%s\n' "$OUTPUT_TREE" > "$SUCCESS_FILE"
else
  rm -f "$LOCK_FILE" "$SUCCESS_FILE"
fi

mirror_results_leaf

echo "ASTRAL-Pro3 finished in ${RUNNING_TIME}s (exit code ${EXIT_CODE}); RF rate: ${RF_RATE}"
echo "Wrote stats to $STAT_FILE"

if [[ "$NO_NOTIFY" == false ]] && command -v curl >/dev/null 2>&1; then
  STATUS_EMOJI=$(if [[ $EXIT_CODE -eq 0 ]]; then echo "✅"; else echo "❌"; fi)
  STATUS_TEXT=$(if [[ $EXIT_CODE -eq 0 ]]; then echo "completed"; else echo "failed (exit $EXIT_CODE)"; fi)
  DISPLAY_OPTS="${ASTRAL_PRO3_OPTS:-<default>}"
  CSV_HEADER="alg,setting,num-taxa,gene-trees,replicate,sb,spmin,spmax,rf-rate,optimal-triplet-score,running-time-s,max-cpu-mb,max-gpu-mb"
  CSV_ROW="astral-pro3,${SETTING_NAME},${TAXA_NUM},${GENE_TREES},${REPLICATE},${SB},${SPMIN},${SPMAX},${RF_RATE},NA,${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB}"
  NOTIFY_BODY="${STATUS_EMOJI} ASTRAL-Pro3 ${STATUS_TEXT} for ${TAXA_NUM} taxa and ${GENE_TREES} gene trees

Setting: ${SETTING_NAME}
Options: ${DISPLAY_OPTS}
RF Rate: ${RF_RATE}
Running time: ${RUNNING_TIME}s
Max CPU RAM: ${MAX_CPU_MB} MB
Max GPU VRAM: ${MAX_GPU_MB} MB

${CSV_HEADER}
${CSV_ROW}

Stats: ${STAT_FILE}"
  curl -s -d "$NOTIFY_BODY" "https://ntfy.sh/${NTFY_CHANNEL_NAME}" >/dev/null 2>&1 || true
fi

exit "$EXIT_CODE"
