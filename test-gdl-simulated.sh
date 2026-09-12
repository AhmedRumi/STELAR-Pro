#!/usr/bin/env bash
# Run one method on an existing GDL replicate. Never regenerate source data.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT}/scripts/phylogeny-data-dir.sh"
source "${ROOT}/scripts/simulated-outputs-mirror.sh"
source "${ROOT}/experiment-setting-name.sh"
INVOCATION=("$0" "$@")
DATA="" DATASET="" REPLICATE=R1 METHOD=stelar-pro OPTS="" BIN="" OUTPUTS=""
FRESH=false NOTIFY=true MIRROR=true TIME=true GPU=true DRY=false

usage() {
  cat <<'EOF'
Usage: ./test-gdl-simulated.sh --dataset taxa10_gt10_dup1_loss0_pop50000000 [options]
  --gdl-data-dir DIR       Default: $PHYLOGENY_DATA_DIR/gdl-simulation/data
  --replicate R1          Existing replicate (default: R1)
  --method METHOD         stelar-pro (default) or astral-pro3
  --opts "..."            Method options, e.g. --cpu or --thread 8
  --astral-pro3-bin FILE  ASTRAL-Pro3 executable override
  --simulated-outputs-dir DIR  Default: $PHYLOGENY_DATA_DIR/outputs/gdl-simulation
  --fresh                 Rerun; preserve previous successful results
  --no-notify             Disable completion notifications
  --no-outputs-mirror     Disable safety copy
  --no-time-monitor      Disable time/RAM monitoring
  --no-gpu-monitor       Disable GPU monitor (always disabled for ASTRAL-Pro3)
  --dry-run               Show selected files and setting without writes
EOF
}
need_value() { [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; exit 2; }; }
while (($#)); do
  case "$1" in
    --dataset) need_value "$@"; DATASET="$2"; shift 2 ;;
    --gdl-data-dir) need_value "$@"; DATA="$2"; shift 2 ;;
    --replicate|-r) need_value "$@"; REPLICATE="$2"; shift 2 ;;
    --method) need_value "$@"; METHOD="$2"; shift 2 ;;
    --opts) need_value "$@"; OPTS="$2"; shift 2 ;;
    --opts=*) OPTS="${1#*=}"; shift ;;
    --astral-pro3-bin) need_value "$@"; BIN="$2"; shift 2 ;;
    --simulated-outputs-dir|--gdl-simulation-outputs-dir) need_value "$@"; OUTPUTS="$2"; shift 2 ;;
    --fresh) FRESH=true; shift ;;
    --no-notify|-nn) NOTIFY=false; shift ;;
    --no-outputs-mirror) MIRROR=false; shift ;;
    --no-time-monitor) TIME=false; shift ;;
    --no-gpu-monitor) GPU=false; shift ;;
    --dry-run) DRY=true; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Error: unknown option: $1" >&2; exit 2 ;;
  esac
done
stelar_pro_gdl_dataset_name_is_valid "$DATASET" || { echo "Error: invalid GDL dataset name: $DATASET" >&2; exit 2; }
[[ "$REPLICATE" =~ ^R[1-9][0-9]*$ ]] || { echo 'Error: replicate must be R1, R2, ...' >&2; exit 2; }
case "$METHOD" in
  stelar-pro) SETTING="$(build_setting_name_from_opts "$OPTS")" ;;
  astral-pro3)
    [[ ! " $OPTS " =~ [[:space:]]--search-space([=[:space:]]) ]] || { echo 'Error: --search-space is STELAR-Pro-only.' >&2; exit 2; }
    SETTING="$(build_astral_pro3_setting_name_from_opts "$OPTS")"; GPU=false ;;
  *) echo "Error: unsupported method: $METHOD" >&2; exit 2 ;;
esac
# Input/output/taxa mapping must stay under this runner's control.
read -r -a TOKENS <<< "$OPTS"
for token in "${TOKENS[@]}"; do
  case "${token%%=*}" in
    -i|--input|-o|--output|-a|--mapping|--gene-species-map|--taxa-file|--extract-taxa|--tag-only|-T|--score|--score-species-tree|--species-tree|-c|--constraint|--log-file)
      echo "Error: $token is not an inference option supported by the GDL runner." >&2; exit 2 ;;
  esac
done
[[ "$SETTING" =~ ^[a-zA-Z0-9][a-zA-Z0-9_.+-]*$ ]] || { echo 'Error: unsafe setting name.' >&2; exit 2; }
DATA="$(stelar_pro_resolve_gdl_data_dir "$DATA")"
[[ "$MIRROR" == false ]] || OUTPUTS="$(stelar_pro_simulated_outputs_dir "$DATA" "$OUTPUTS")"
RUN="${DATA}/${DATASET}/${REPLICATE}"
INPUT="${RUN}/true-genetrees/all_gt.trees"
REFERENCE="${RUN}/species-tree/s_tree.trees"
for file in "$INPUT" "$REFERENCE"; do
  [[ -s "$file" ]] || { echo "Error: required combined input/reference missing or empty: $file" >&2; exit 6; }
done
PARENT="${RUN}/${METHOD}-outputs"
RESULT="${PARENT}/${SETTING}"
[[ "$(realpath -m -- "$PARENT")" == "$PARENT" && ! -L "$RESULT" ]] || {
  echo 'Error: refusing symlinked dataset/replicate/result directories.' >&2; exit 2;
}
printf '%s %s/%s [%s]\n  Input: %s\n  Reference: %s\n  Results: %s\n' "$METHOD" "$DATASET" "$REPLICATE" "$SETTING" "$INPUT" "$REFERENCE" "$RESULT"
[[ "$DRY" == false ]] || exit 0
mkdir -p -- "$PARENT"
[[ ! -L "${PARENT}/.${SETTING}.running.lock" ]] || { echo 'Error: refusing symlinked running lock.' >&2; exit 2; }
exec 9>>"${PARENT}/.${SETTING}.running.lock"
flock -n 9 || { echo "Error: this replicate/method/setting is already running." >&2; exit 4; }
mirror_leaf() {
  [[ "$MIRROR" == false ]] || stelar_pro_mirror_simulated_results "$DATA" "$OUTPUTS" "$1"
}
if [[ "$FRESH" == false && -s "${RESULT}/out-${METHOD}.tre" && -f "${RESULT}/.${METHOD}.success" && -s "${RESULT}/stat-${METHOD}.csv" ]]; then
  if python3 - "${RESULT}/input-preparation.json" "$INPUT" "$REFERENCE" <<'PY'
import hashlib, json, sys
def digest(path):
    value = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()
try:
    with open(sys.argv[1]) as stream:
        metadata = json.load(stream)
    valid = metadata['input_sha256'] == digest(sys.argv[2]) and metadata['reference_sha256'] == digest(sys.argv[3])
except (OSError, ValueError, KeyError):
    valid = False
sys.exit(0 if valid else 1)
PY
  then
    mirror_leaf "$RESULT"
    echo 'SKIPPING: completed setting and input fingerprints match; use --fresh to rerun.'
    exit 0
  fi
  echo 'Input fingerprints changed/missing: rerunning and preserving the previous result.'
fi

# Normalize into a temporary directory, never inside a mirrored result leaf.
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelar-gdl-input.XXXXXX")"
STAGE="$(mktemp -d "${PARENT}/.${SETTING}.run.XXXXXX")"
trap 'rm -rf -- "$WORK"' EXIT
NORMALIZED="${WORK}/species-labelled.trees"
PREPARE=(python3 "${ROOT}/scripts/prepare_gdl_input.py" --input "$INPUT" --reference "$REFERENCE"
  --output "$NORMALIZED" --metadata "${STAGE}/input-preparation.json")
IFS=_ read -r TAXA GT DUP LOSS POP <<< "$DATASET"
TAXA="${TAXA#taxa}"; GT="${GT#gt}"; DUP="${DUP#dup}"; LOSS="${LOSS#loss}"; POP="${POP#pop}"
PREPARE+=(--expected-genes "$GT")
TREE="${STAGE}/out-${METHOD}.tre"
if [[ "$METHOD" == stelar-pro ]]; then
  CMD=(bash "${ROOT}/run-stelar-pro-with-monitor.sh" -i "$NORMALIZED" -o "$TREE" --stelar-pro-root "$ROOT" --no-notify)
else
  CMD=(bash "${ROOT}/run-astral-pro3-with-monitor.sh" --input "$NORMALIZED" --output "$TREE" --no-notify)
  [[ -z "$BIN" ]] || CMD+=(--astral-pro3-bin "$BIN")
fi
[[ -z "$OPTS" ]] || CMD+=(--opts "$OPTS")
[[ "$TIME" == true ]] || CMD+=(--no-time-monitor)
[[ "$GPU" == true ]] || CMD+=(--no-gpu-monitor)
EXIT_CODE=0
if "${PREPARE[@]}" >"${STAGE}/input-preparation.log" 2>&1; then
  cat "${STAGE}/input-preparation.log"
  set +e
  "${CMD[@]}" 2>&1 | tee "${STAGE}/.${METHOD}_run.log"
  PIPE_EXIT=("${PIPESTATUS[@]}")
  EXIT_CODE=${PIPE_EXIT[0]}
  set -e
  [[ "$EXIT_CODE" != 0 || "${PIPE_EXIT[1]}" == 0 ]] || EXIT_CODE=7
  [[ "$EXIT_CODE" != 0 || -s "$TREE" ]] || EXIT_CODE=5
else
  EXIT_CODE=6
  cat "${STAGE}/input-preparation.log" >&2
fi
RUNTIME=NA CPU=NA GPU_MB=NA SCORE=NA RF=NA RF_STATUS=unavailable
WRAPPER_STATS="${TREE%.tre}_stats.csv"
if [[ -s "$WRAPPER_STATS" ]]; then
  IFS=, read -r _ _ _ RUNTIME CPU GPU_MB SCORE _ _ < <(sed -n '2p' "$WRAPPER_STATS")
fi
if [[ "$EXIT_CODE" == 0 ]]; then
  if python3 "${ROOT}/scripts/calculate_rf_rate.py" "$TREE" "$REFERENCE" >"${STAGE}/rf.log" 2>&1; then
    RF="$(awk '/^RF rate:/ {print $3}' "${STAGE}/rf.log")"; RF_STATUS=ok
  else
    RF_STATUS=error
    echo 'Warning: strict RF unavailable; see rf.log (no pruning or label mismatch is silently accepted).' >&2
  fi
fi
HEADER='alg,setting,num-taxa,gene-trees,replicate,dup,loss,pop,rf-rate,optimal-triplet-score,running-time-s,max-cpu-mb,max-gpu-mb,dataset,exit-code,rf-status'
ROW="${METHOD},${SETTING},${TAXA},${GT},${REPLICATE},${DUP},${LOSS},${POP},${RF},${SCORE},${RUNTIME},${CPU},${GPU_MB},${DATASET},${EXIT_CODE},${RF_STATUS}"
printf '%s\n%s\n' "$HEADER" "$ROW" >"${STAGE}/stat-${METHOD}.csv"
COMMAND="${TREE%.tre}.command"
{
  echo
  echo '# GDL run context (normalization is required before replaying inference)'
  echo "# original input: $INPUT"
  echo "# reference species tree: $REFERENCE"
  echo "# dataset: $DATASET; replicate: $REPLICATE; setting: $SETTING"
  echo "# revision: $(stelar_pro_git_revision "$ROOT")"
  echo '# input-preparation.json contains source/reference/normalized SHA-256 IDs'
  printf '# invocation: '; stelar_pro_print_shell_command "${INVOCATION[@]}"
  printf '# preparation command: '; stelar_pro_print_shell_command "${PREPARE[@]}"
  printf '# monitor command: '; stelar_pro_print_shell_command "${CMD[@]}"
  echo "# exit code: $EXIT_CODE; RF rate: $RF; RF status: $RF_STATUS"
} >>"$COMMAND"
STAMP="$(date -u +%Y%m%dT%H%M%S).$$"
if [[ "$EXIT_CODE" == 0 ]]; then
  if [[ -e "$RESULT" ]]; then
    # Keep the old leaf as an explicit archive, not a misleading current result.
    ARCHIVE="${PARENT}/${SETTING}__previous_${STAMP}"
    mv -- "$RESULT" "$ARCHIVE"
    if ! mv -- "$STAGE" "$RESULT"; then mv -- "$ARCHIVE" "$RESULT"; exit 1; fi
  else
    mv -- "$STAGE" "$RESULT"
  fi
  FINAL="$RESULT"
else
  FINAL="${PARENT}/${SETTING}__failed_${STAMP}"
  mv -- "$STAGE" "$FINAL"
fi
echo "# final result directory: $FINAL" >>"${FINAL}/out-${METHOD}.command"
# The executed monitor command intentionally retains its exact temporary paths.
# Give users a directly runnable recipe that recreates normalization safely.
REPLAY=(bash "${ROOT}/test-gdl-simulated.sh" --gdl-data-dir "$DATA" --dataset "$DATASET"
  --replicate "$REPLICATE" --method "$METHOD" --opts "$OPTS" --fresh --no-notify)
[[ -z "$BIN" ]] || REPLAY+=(--astral-pro3-bin "$BIN")
[[ "$TIME" == true ]] || REPLAY+=(--no-time-monitor)
[[ "$GPU" == true ]] || REPLAY+=(--no-gpu-monitor)
if [[ "$MIRROR" == true ]]; then REPLAY+=(--simulated-outputs-dir "$OUTPUTS"); else REPLAY+=(--no-outputs-mirror); fi
{
  printf '#!/usr/bin/env bash\nset -euo pipefail\n'
  stelar_pro_print_shell_command "${REPLAY[@]}"
} >"${FINAL}/replay-gdl.sh"
# Sidecars describe persistent source/result paths after the staging directory
# has been renamed; the command record still describes the exact execution.
if [[ -f "${FINAL}/out-${METHOD}_stats.csv" ]]; then
  python3 - "${FINAL}/out-${METHOD}_stats.csv" "$INPUT" "${FINAL}/out-${METHOD}.tre" "$RF" <<'PY'
import csv, sys
with open(sys.argv[1], newline='') as stream:
    rows = list(csv.reader(stream))
if len(rows) >= 2:
    columns = {name: i for i, name in enumerate(rows[0])}
    for name, value in (('input_file', sys.argv[2]), ('output_file', sys.argv[3]), ('rf_rate', sys.argv[4])):
        rows[1][columns[name]] = value
    with open(sys.argv[1], 'w', newline='') as stream:
        csv.writer(stream).writerows(rows)
PY
fi
if [[ "$EXIT_CODE" == 0 ]]; then
  touch "${FINAL}/.${METHOD}.lock"
  printf 'exit_code=0\noutput=%s\n' "${FINAL}/out-${METHOD}.tre" >"${FINAL}/.${METHOD}.success"
fi
[[ -z "${ARCHIVE:-}" ]] || mirror_leaf "$ARCHIVE"
mirror_leaf "$FINAL"
echo "Finished ${METHOD}: exit ${EXIT_CODE}; RF ${RF}; time ${RUNTIME}s; stats: ${FINAL}/stat-${METHOD}.csv"
if [[ "$NOTIFY" == true ]] && command -v curl >/dev/null 2>&1; then
  STATUS='✅ completed'; [[ "$EXIT_CODE" == 0 ]] || STATUS="❌ failed (exit ${EXIT_CODE})"
  BODY="${STATUS}: ${METHOD} on GDL simulated data

Dataset: ${DATASET}
Replicate: ${REPLICATE}
Taxa: ${TAXA}; gene trees: ${GT}; dup: ${DUP}; loss: ${LOSS}; pop: ${POP}
Setting: ${SETTING}
Options: ${OPTS:-<default>}
RF Rate: ${RF} (${RF_STATUS})
Optimal triplet score: ${SCORE}
Running time: ${RUNTIME}s
Max CPU RAM: ${CPU} MB
Max GPU VRAM: ${GPU_MB} MB

${HEADER}
${ROW}

Stats: ${FINAL}/stat-${METHOD}.csv"
  curl --max-time 15 -s -d "$BODY" "https://ntfy.sh/${NTFY_CHANNEL_NAME:-anik-phylo-stx}" >/dev/null 2>&1 || true
fi
exit "$EXIT_CODE"
