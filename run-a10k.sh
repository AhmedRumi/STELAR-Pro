#!/usr/bin/env bash
# run-a10k.sh
# STELAR-Pro runner for the 10k-simphy dataset layout.

set -euo pipefail

NTFY_CHANNEL_NAME="${NTFY_CHANNEL_NAME:-anik-phylo-stx}"

TREE_TYPES_RAW="estimated"
DATA_DIR=""
REPLICATES_SPEC=""
START_REP=""
END_REP=""
FRESH=false
STELAR_PRO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STELAR_PRO_OPTS="-vv"
STELAR_PRO_OPTS_LIST_RAW=""
TIME_MONITOR=true
GPU_MONITOR=true
NO_NOTIFY=false

source "${STELAR_PRO_ROOT}/experiment-setting-name.sh"

csv_get_field() {
  local file="$1"
  shift
  local header data
  header="$(head -n1 "$file" 2>/dev/null || true)"
  data="$(sed -n '2p' "$file" 2>/dev/null || true)"
  if [[ -z "$header" || -z "$data" ]]; then
    echo ""
    return 0
  fi
  IFS=',' read -r -a headers <<< "$header"
  IFS=',' read -r -a values <<< "$data"
  for key in "$@"; do
    for i in "${!headers[@]}"; do
      if [[ "${headers[$i]}" == "$key" ]]; then
        echo "${values[$i]:-}"
        return 0
      fi
    done
  done
  echo ""
}

# Verbosity flags such as -v/-vv are ignored when constructing setting names.

print_help() {
  cat <<EOF
run-a10k.sh

Required:
  --data-dir           Path to A10K dataset root containing 10k-simphy/

Optional:
  --tree-type          estimated | true, or a semicolon-separated list
                       such as "true;estimated" (default: estimated)
  --replicates         Replicates to run, e.g. "1-20" or "R1,R2"
  --start-rep, -sr     Start replicate number
  --end-rep, -er       End replicate number
  --stelar-pro-root       Path to STELAR-Pro root
  --opts, --alg-opts   Extra options for the selected algorithm
  --opts-list, --alg-opts-list
                       Semicolon-separated list of option strings to loop over
  --fresh              Force rerun even if stat-stelar-pro.csv exists
  --no-time-monitor    Disable time monitoring
  --no-gpu-monitor     Disable GPU monitoring
  --no-notify, -nn     Disable ntfy notifications

Examples:
  ./run-a10k.sh --data-dir /path/to/10k-astral-dataset --tree-type estimated --opts "--threads 8 -vv"
  ./run-a10k.sh --data-dir /path/to/10k-astral-dataset --tree-type "true;estimated" --opts "--gpu-strict -vv"
  The first example setting is threads_8.
  Verbosity is ignored; other meaningful options are appended to the name.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA_DIR="$2"; shift 2 ;;
    --tree-type) TREE_TYPES_RAW="$2"; shift 2 ;;
    --tree-type=*) TREE_TYPES_RAW="${1#*=}"; shift ;;
    --replicates) REPLICATES_SPEC="$2"; shift 2 ;;
    --start-rep|-sr) START_REP="$2"; shift 2 ;;
    --end-rep|-er) END_REP="$2"; shift 2 ;;
    --stelar-pro-root|--stelar-root) STELAR_PRO_ROOT="$2"; shift 2 ;;
    --opts|--alg-opts|--stelar-pro-opts|--stelar-opts) STELAR_PRO_OPTS="$2"; shift 2 ;;
    --opts-list|--alg-opts-list|--stelar-pro-opts-list|--stelar-opts-list) STELAR_PRO_OPTS_LIST_RAW="$2"; shift 2 ;;
    --fresh) FRESH=true; shift ;;
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --help|-h) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [[ -z "$DATA_DIR" ]]; then
  echo "Error: --data-dir is required."
  exit 2
fi

STELAR_PRO_ROOT="$(realpath "$STELAR_PRO_ROOT")"
PYTHON_BIN="${STELAR_PRO_PYTHON:-${STELAR_PRO_ROOT}/.venv/bin/python}"
[[ -x "$PYTHON_BIN" ]] || PYTHON_BIN="python3"
DATA_DIR="$(realpath "$DATA_DIR")"
SIMPHY_DIR="${DATA_DIR%/}/10k-simphy"
if [[ ! -d "$SIMPHY_DIR" ]]; then
  echo "Error: expected 10k-simphy at $SIMPHY_DIR"
  exit 3
fi

TREE_TYPES=()
IFS=';' read -r -a raw_tree_types <<< "$TREE_TYPES_RAW"
for tree_type in "${raw_tree_types[@]}"; do
  tree_type="${tree_type//[[:space:]]/}"
  tree_type="${tree_type,,}"
  [[ -n "$tree_type" ]] || continue
  case "$tree_type" in
    true|estimated) ;;
    *)
      echo "Error: invalid --tree-type value '$tree_type' (expected true, estimated, or a semicolon-separated list)."
      exit 2
      ;;
  esac

  duplicate=false
  for existing_tree_type in "${TREE_TYPES[@]}"; do
    if [[ "$existing_tree_type" == "$tree_type" ]]; then
      duplicate=true
      break
    fi
  done
  [[ "$duplicate" == false ]] && TREE_TYPES+=("$tree_type")
done
if [[ ${#TREE_TYPES[@]} -eq 0 ]]; then
  echo "Error: --tree-type must contain at least one of: true, estimated."
  exit 2
fi

STELAR_PRO_OPTS_LIST=()
if [[ -n "$STELAR_PRO_OPTS_LIST_RAW" ]]; then
  IFS=';' read -r -a raw_opts_list <<< "$STELAR_PRO_OPTS_LIST_RAW"
  for opts in "${raw_opts_list[@]}"; do
    opts="$(echo "$opts" | sed 's/^ *//;s/ *$//')"
    [[ -n "$opts" ]] && STELAR_PRO_OPTS_LIST+=("$opts")
  done
fi
if [[ ${#STELAR_PRO_OPTS_LIST[@]} -eq 0 ]]; then
  STELAR_PRO_OPTS_LIST+=("${STELAR_PRO_OPTS}")
fi
echo "[DEBUG] opts list (${#STELAR_PRO_OPTS_LIST[@]} items): ${STELAR_PRO_OPTS_LIST[*]}"
echo "[DEBUG] tree types (${#TREE_TYPES[@]} items): ${TREE_TYPES[*]}"
echo "[DEBUG] replicates spec: '${REPLICATES_SPEC}' | fresh: ${FRESH}"

REPL_LIST=()
if [[ -n "$START_REP" || -n "$END_REP" ]]; then
  for i in $(seq "$START_REP" "$END_REP"); do REPL_LIST+=("R${i}"); done
elif [[ -n "$REPLICATES_SPEC" ]]; then
  if [[ "$REPLICATES_SPEC" =~ ^[0-9]+-[0-9]+$ ]]; then
    start="${REPLICATES_SPEC%-*}"
    end="${REPLICATES_SPEC#*-}"
    for i in $(seq "$start" "$end"); do REPL_LIST+=("R${i}"); done
  else
    IFS=',' read -r -a parts <<< "$REPLICATES_SPEC"
    for p in "${parts[@]}"; do
      p="${p// /}"
      [[ "$p" =~ ^R ]] || p="R${p}"
      REPL_LIST+=("$p")
    done
  fi
else
  while IFS= read -r -d '' d; do
    REPL_LIST+=("$(basename "$d")")
  done < <(find "$SIMPHY_DIR" -maxdepth 1 -type d -name 'R*' -print0 | sort -z -V)
fi

echo "[DEBUG] replicate list (${#REPL_LIST[@]} items): ${REPL_LIST[*]}"

for TREE_TYPE in "${TREE_TYPES[@]}"; do
  echo "==> Processing A10K tree type: ${TREE_TYPE}"
  for REPL in "${REPL_LIST[@]}"; do
  REPL_DIR="${SIMPHY_DIR%/}/${REPL}"
  if [[ ! -d "$REPL_DIR" ]]; then
    echo "[DEBUG] SKIP ${REPL}: directory not found: ${REPL_DIR}"
    continue
  fi

  if [[ "$TREE_TYPE" == "estimated" ]]; then
    GT_DIR="${REPL_DIR}/estimatedgenetrees"
    GT_FILE="${GT_DIR}/estimatedgenetrees.tre"
    ROOTED_GT="${GT_DIR}/estimatedgenetrees.rooted.tre"
    if [[ ! -f "$ROOTED_GT" ]]; then
      if [[ ! -x "${STELAR_PRO_ROOT%/}/process_unrooted.sh" ]]; then
        echo "Error: process_unrooted.sh not found or not executable at ${STELAR_PRO_ROOT%/}/process_unrooted.sh"
        exit 7
      fi
      echo "Rooting estimated gene trees for ${REPL} with outgroup 0..."
      "${STELAR_PRO_ROOT%/}/process_unrooted.sh" -i "$GT_FILE" -o "$ROOTED_GT" -og "0"
    fi
    GT_FILE="$ROOTED_GT"
  else
    GT_FILE="${REPL_DIR}/truegenetrees"
  fi
  TRUE_TREE="${REPL_DIR}/s_tree.trees"
  if [[ ! -f "$GT_FILE" || ! -f "$TRUE_TREE" ]]; then
    echo "[DEBUG] SKIP ${REPL}: missing files (gt=${GT_FILE} exists=$([ -f "$GT_FILE" ] && echo yes || echo no), true_tree=${TRUE_TREE} exists=$([ -f "$TRUE_TREE" ] && echo yes || echo no))"
    continue
  fi

  for STELAR_PRO_OPTS_ITEM in "${STELAR_PRO_OPTS_LIST[@]}"; do
    SETTING_NAME="$(build_setting_name_from_opts "$STELAR_PRO_OPTS_ITEM")"
    OUT_DIR="${REPL_DIR}/stelar-pro-outputs/${TREE_TYPE}/${SETTING_NAME}"
    OUT_FILE="${OUT_DIR}/out-stelar-pro.tre"
    STAT_FILE="${OUT_DIR}/stat-stelar-pro.csv"

    if [[ "$FRESH" == false && -f "$STAT_FILE" ]]; then
      echo "SKIPPING: ${STAT_FILE} exists."
      continue
    elif [[ "$FRESH" == true && -f "$STAT_FILE" ]]; then
      echo "[DEBUG] --fresh set, overwriting existing: ${STAT_FILE}"
    fi

    mkdir -p "$OUT_DIR"
    CMD=("${STELAR_PRO_ROOT}/run-stelar-pro-with-monitor.sh" -i "$GT_FILE" -o "$OUT_FILE" --stelar-pro-root "$STELAR_PRO_ROOT")
    if [[ "$TIME_MONITOR" == false ]]; then CMD+=(--no-time-monitor); fi
    if [[ "$GPU_MONITOR" == false ]]; then CMD+=(--no-gpu-monitor); fi
    if [[ "$NO_NOTIFY" == true ]]; then CMD+=(--no-notify); fi
    if [[ -n "$STELAR_PRO_OPTS_ITEM" ]]; then
      CMD+=(--opts "$STELAR_PRO_OPTS_ITEM")
    fi

    echo "==> Running stelar-pro on ${REPL} (${TREE_TYPE}, ${SETTING_NAME})"
    echo "Command: ${CMD[*]}"
    set +e
    "${CMD[@]}"
    RUN_EXIT=$?
    set -e

    SIDE_STATS="${OUT_FILE%.tre}_stats.csv"
    if [[ "$RUN_EXIT" -ne 0 || ! -f "$SIDE_STATS" ]]; then
      echo "Run failed for ${REPL} (${TREE_TYPE}, ${SETTING_NAME}); skipping RF/stat summary."
      continue
    fi

    RUNNING_TIME="$(csv_get_field "$SIDE_STATS" "running_time_s" "running-time-s")"
    MAX_CPU_MB="$(csv_get_field "$SIDE_STATS" "max_cpu_mb" "max-cpu-mb")"
    MAX_GPU_MB="$(csv_get_field "$SIDE_STATS" "max_gpu_mb" "max-gpu-mb")"
    OPTIMAL_TRIPLET_SCORE="$(csv_get_field "$SIDE_STATS" "optimal_triplet_score" "optimal-triplet-score")"
    EXIT_CODE="$(csv_get_field "$SIDE_STATS" "exit_code" "exit-code")"
    if [[ -z "$EXIT_CODE" ]]; then
      EXIT_CODE="$RUN_EXIT"
    fi

    RF_RATE="NA"
    if [[ -f "$OUT_FILE" && -f "$TRUE_TREE" ]]; then
      rf_output=$("$PYTHON_BIN" "${STELAR_PRO_ROOT}/rf.py" "$OUT_FILE" "$TRUE_TREE" 2>&1) || true
      rf_line=$(echo "$rf_output" | grep -i "Robinson-Foulds distance" | tail -n1 || true)
      if [[ -n "$rf_line" ]]; then
        RF_RATE=$(echo "$rf_line" | grep -Eo '[0-9]+(\.[0-9]+)?' | tail -n1 || echo "NA")
      fi
    fi

    echo "alg,setting,replicate,tree_type,rf-rate,optimal-triplet-score,running-time-s,max-cpu-mb,max-gpu-mb" > "$STAT_FILE"
    echo "stelar-pro,${SETTING_NAME},${REPL},${TREE_TYPE},${RF_RATE},${OPTIMAL_TRIPLET_SCORE},${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB}" >> "$STAT_FILE"
    echo
    echo "=== A10K STELAR-Pro Summary ==="
    echo "Replicate:      ${REPL}"
    echo "Tree type:      ${TREE_TYPE}"
    echo "Setting:        ${SETTING_NAME}"
    echo "RF rate:        ${RF_RATE}"
    echo "Triplet score:  ${OPTIMAL_TRIPLET_SCORE}"
    echo "Running time:   ${RUNNING_TIME}s"
    echo "Max CPU RAM:    ${MAX_CPU_MB} MB"
    echo "Max GPU VRAM:   ${MAX_GPU_MB} MB"
    echo "Output tree:    ${OUT_FILE}"
    echo "Stats file:     ${STAT_FILE}"
    echo "Saved $STAT_FILE"

    if [[ "$NO_NOTIFY" == false ]] && command -v curl >/dev/null 2>&1; then
      curl -s -d "✅ STELAR-Pro A10K completed

Replicate: ${REPL}
Tree type: ${TREE_TYPE}
Setting: ${SETTING_NAME}

RF: ${RF_RATE}
Triplet score: ${OPTIMAL_TRIPLET_SCORE}
Time: ${RUNNING_TIME}s
CPU: ${MAX_CPU_MB} MB
GPU: ${MAX_GPU_MB} MB
Exit: ${EXIT_CODE}

Tree: $(basename "$OUT_FILE")
Stats: $(basename "$STAT_FILE")" "https://ntfy.sh/${NTFY_CHANNEL_NAME}" >/dev/null 2>&1 || true
    fi
    done
  done
done
