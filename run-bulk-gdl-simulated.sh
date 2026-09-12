#!/usr/bin/env bash
# Sweep EXISTING actual gene-duplication/loss simulations (no generation).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT}/scripts/phylogeny-data-dir.sh"
source "${ROOT}/scripts/simulated-outputs-mirror.sh"
source "${ROOT}/experiment-setting-name.sh"

# Edit these lists, or override them on the command line (comma/space separated).
TAXA_LIST=(10)
GT_LIST=(10)
DUP_LIST=(1)
LOSS_LIST=(0)
POP_LIST=(50000000)
NUM_REPLICATES=1
METHODS_RAW=stelar-pro
OPTS_LIST=("")
STELAR_OPTS_RAW="" ASTRAL_OPTS_RAW="" STELAR_OPTS_SET=false ASTRAL_OPTS_SET=false
DATA="" OUTPUTS="" BIN=""
DRY=false YES=false STRICT=false
FORWARDED=()
usage() {
  cat <<'EOF'
Usage: ./run-bulk-gdl-simulated.sh [options]

Run existing <dataset>/R*/true-genetrees/all_gt.trees against
<dataset>/R*/species-tree/s_tree.trees. Separate gene trees are never required.

  --taxa-list LIST        Default: 10
  --gt-list LIST          Default: 10 (alias: --genes-list)
  --dup-list LIST         Default: 1 (directory's dup parameter, not a rate)
  --loss-list LIST        Default: 0
  --pop-list LIST         Default: 50000000
  --num-replicates N      Select R1 through RN (default: 1)
  --method METHOD        stelar-pro (default) or astral-pro3
  --methods LIST         Run both: "stelar-pro,astral-pro3"
  --opts "..."           One method option string
  --opts-list LIST        Semicolon-separated option strings
  --stelar-pro-opts-list LIST   Method-specific sweep options
  --astral-pro3-opts-list LIST  Method-specific sweep options
  --gdl-data-dir DIR     Default: $PHYLOGENY_DATA_DIR/gdl-simulation/data
  --simulated-outputs-dir DIR  Default: $PHYLOGENY_DATA_DIR/outputs/gdl-simulation
  --astral-pro3-bin FILE Executable override
  --strict-missing        Abort before inference if any selected input is missing
  --fresh                Rerun completed settings, preserving old results
  --no-notify            Disable notifications (otherwise enabled)
  --no-outputs-mirror    Disable safety copy
  --no-time-monitor     Disable time/RAM monitoring
  --no-gpu-monitor      Disable STELAR GPU monitoring; ASTRAL is CPU-only
  --dry-run              Print plan only, without writes or inference
  --yes, -y              Skip interactive confirmation

Lists accept commas or spaces; directory parameter spelling is matched exactly.
Missing inputs are skipped with warnings by default. No simulations are created.
EOF
}
need_value() { [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; exit 2; }; }
parse_list() {
  local raw="${2//,/ }"
  read -r -a "$1" <<< "$raw"
}
parse_opts() {
  # Preserve an empty string as the default method setting.
  local -n destination="$1"
  if [[ -z "$2" ]]; then destination=(""); else IFS=';' read -r -a destination <<< "$2"; fi
}
while (($#)); do
  case "$1" in
    --taxa-list) need_value "$@"; parse_list TAXA_LIST "$2"; shift 2 ;;
    --gt-list|--genes-list|--gene-trees-list) need_value "$@"; parse_list GT_LIST "$2"; shift 2 ;;
    --dup-list) need_value "$@"; parse_list DUP_LIST "$2"; shift 2 ;;
    --loss-list) need_value "$@"; parse_list LOSS_LIST "$2"; shift 2 ;;
    --pop-list) need_value "$@"; parse_list POP_LIST "$2"; shift 2 ;;
    --num-replicates) need_value "$@"; NUM_REPLICATES="$2"; shift 2 ;;
    --method|--methods) need_value "$@"; METHODS_RAW="$2"; shift 2 ;;
    --opts) need_value "$@"; OPTS_LIST=("$2"); shift 2 ;;
    --opts=*) OPTS_LIST=("${1#*=}"); shift ;;
    --opts-list) need_value "$@"; parse_opts OPTS_LIST "$2"; shift 2 ;;
    --opts-list=*) parse_opts OPTS_LIST "${1#*=}"; shift ;;
    --stelar-pro-opts-list) need_value "$@"; STELAR_OPTS_RAW="$2"; STELAR_OPTS_SET=true; shift 2 ;;
    --astral-pro3-opts-list) need_value "$@"; ASTRAL_OPTS_RAW="$2"; ASTRAL_OPTS_SET=true; shift 2 ;;
    --gdl-data-dir) need_value "$@"; DATA="$2"; shift 2 ;;
    --simulated-outputs-dir|--gdl-simulation-outputs-dir) need_value "$@"; OUTPUTS="$2"; shift 2 ;;
    --astral-pro3-bin) need_value "$@"; BIN="$2"; shift 2 ;;
    --strict-missing) STRICT=true; shift ;;
    --fresh|--no-notify|-nn|--no-outputs-mirror|--no-time-monitor|--no-gpu-monitor) FORWARDED+=("$1"); shift ;;
    --dry-run) DRY=true; shift ;;
    --yes|-y) YES=true; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Error: unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ "$NUM_REPLICATES" =~ ^[1-9][0-9]*$ ]] || { echo 'Error: --num-replicates must be positive.' >&2; exit 2; }
for name in TAXA_LIST GT_LIST POP_LIST DUP_LIST LOSS_LIST; do
  declare -n values="$name"
  ((${#values[@]})) || { echo "Error: $name is empty." >&2; exit 2; }
  for value in "${values[@]}"; do
    if [[ "$name" == DUP_LIST || "$name" == LOSS_LIST ]]; then
      [[ "$value" =~ ^[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$ ]] || { echo "Error: invalid $name value: $value" >&2; exit 2; }
    else
      [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "Error: invalid $name value: $value" >&2; exit 2; }
    fi
  done
  unset -n values
done
parse_list METHODS "$METHODS_RAW"
((${#METHODS[@]})) || { echo 'Error: empty method list.' >&2; exit 2; }
STELAR_OPTIONS=("${OPTS_LIST[@]}"); ASTRAL_OPTIONS=("${OPTS_LIST[@]}")
[[ "$STELAR_OPTS_SET" == false ]] || parse_opts STELAR_OPTIONS "$STELAR_OPTS_RAW"
[[ "$ASTRAL_OPTS_SET" == false ]] || parse_opts ASTRAL_OPTIONS "$ASTRAL_OPTS_RAW"
for method in "${METHODS[@]}"; do
  case "$method" in
    stelar-pro) OPTIONS=("${STELAR_OPTIONS[@]}") ;;
    astral-pro3) OPTIONS=("${ASTRAL_OPTIONS[@]}") ;;
    *) echo "Error: unsupported method: $method" >&2; exit 2 ;;
  esac
  for opts in "${OPTIONS[@]}"; do
    if [[ "$method" == astral-pro3 && " $opts " =~ [[:space:]]--search-space([=[:space:]]) ]]; then
      echo 'Error: --search-space is STELAR-only; use method-specific option lists.' >&2; exit 2
    fi
  done
done
DATA="$(stelar_pro_resolve_gdl_data_dir "$DATA")"
COMMANDS_DATASET=() COMMANDS_REPLICATE=() COMMANDS_METHOD=() COMMANDS_OPTS=()
MISSING=0
for taxa in "${TAXA_LIST[@]}"; do
for gt in "${GT_LIST[@]}"; do
for dup in "${DUP_LIST[@]}"; do
for loss in "${LOSS_LIST[@]}"; do
for pop in "${POP_LIST[@]}"; do
  dataset="taxa${taxa}_gt${gt}_dup${dup}_loss${loss}_pop${pop}"
  for ((r=1; r<=NUM_REPLICATES; r++)); do
    replicate="R${r}"; run="${DATA}/${dataset}/${replicate}"
    if [[ ! -s "${run}/true-genetrees/all_gt.trees" || ! -s "${run}/species-tree/s_tree.trees" ]]; then
      echo "Warning: skipping ${dataset}/${replicate}: required combined gene trees or species tree missing/empty." >&2
      MISSING=$((MISSING + 1)); continue
    fi
    for method in "${METHODS[@]}"; do
      if [[ "$method" == stelar-pro ]]; then OPTIONS=("${STELAR_OPTIONS[@]}"); else OPTIONS=("${ASTRAL_OPTIONS[@]}"); fi
      for opts in "${OPTIONS[@]}"; do
        COMMANDS_DATASET+=("$dataset"); COMMANDS_REPLICATE+=("$replicate")
        COMMANDS_METHOD+=("$method"); COMMANDS_OPTS+=("$opts")
      done
    done
  done
done; done; done; done; done
[[ "$STRICT" == false || "$MISSING" == 0 ]] || { echo "Error: ${MISSING} missing selection(s); no methods ran." >&2; exit 6; }
((${#COMMANDS_DATASET[@]})) || { echo 'Error: no existing replicates match the selected lists.' >&2; exit 6; }
build_command() {
  local i="$1"
  CMD=(bash "${ROOT}/test-gdl-simulated.sh" --gdl-data-dir "$DATA" --dataset "${COMMANDS_DATASET[$i]}"
    --replicate "${COMMANDS_REPLICATE[$i]}" --method "${COMMANDS_METHOD[$i]}" --opts "${COMMANDS_OPTS[$i]}" "${FORWARDED[@]}")
  [[ -z "$OUTPUTS" ]] || CMD+=(--simulated-outputs-dir "$OUTPUTS")
  [[ -z "$BIN" ]] || CMD+=(--astral-pro3-bin "$BIN")
}
echo "GDL plan: ${#COMMANDS_DATASET[@]} method/setting run(s); ${MISSING} missing replicate selection(s)."
# Preflight every method/setting before starting expensive inference.
for i in "${!COMMANDS_DATASET[@]}"; do
  build_command "$i"
  "${CMD[@]}" --dry-run
done
[[ "$DRY" == false ]] || { echo 'Dry run complete; no files changed and no inference ran.'; exit 0; }
if [[ "$YES" == false && -t 0 && -t 1 ]]; then
  read -r -p 'Run this GDL plan? [y/N] ' answer
  [[ "$answer" =~ ^[Yy]([Ee][Ss])?$ ]] || { echo 'Cancelled.'; exit 0; }
fi
FAILED=0
for i in "${!COMMANDS_DATASET[@]}"; do
  build_command "$i"
  if ! "${CMD[@]}"; then FAILED=$((FAILED + 1)); fi
done
echo "GDL sweep finished: $((${#COMMANDS_DATASET[@]} - FAILED)) successful/already-completed; ${FAILED} failed; ${MISSING} missing selections skipped."
[[ "$FAILED" == 0 ]]
