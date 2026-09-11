#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelar-pro-astral-pro3-test.XXXXXX")"
trap 'status=$?; rm -rf -- "$WORK"; exit "$status"' EXIT

FAKE_BIN="${WORK}/astral-pro3"
FAKE_LOG="${WORK}/args.log"
export FAKE_LOG

printf '%s\n' \
  '#!/usr/bin/env bash' \
  'set -euo pipefail' \
  'printf "%s\n" "$@" > "$FAKE_LOG"' \
  'input=""; output=""' \
  'while [[ $# -gt 0 ]]; do' \
  '  case "$1" in' \
  '    -i|--input) input="$2"; shift 2 ;;' \
  '    -o|--output) output="$2"; shift 2 ;;' \
  '    --thread|--seed|--verbose) shift 2 ;;' \
  '    *) shift ;;' \
  '  esac' \
  'done' \
  '[[ -f "$input" && -n "$output" ]]' \
  'printf "((a,b),(c,d));\n" > "$output"' \
  > "$FAKE_BIN"
chmod +x "$FAKE_BIN"

# The arbitrary-file wrapper forwards native options and produces its sidecar.
INPUT="${WORK}/genes.tre"
OUTPUT="${WORK}/direct.tre"
printf '((a,b),(c,d));\n' > "$INPUT"
"${ROOT}/run-astral-pro3-with-monitor.sh" \
  -i "$INPUT" -o "$OUTPUT" --astral-pro3-bin "$FAKE_BIN" \
  --opts '--thread 3 --seed 7 --verbose 1' \
  --no-time-monitor --no-gpu-monitor --no-notify >/dev/null
[[ -s "$OUTPUT" && -s "${WORK}/direct_stats.csv" ]]
grep -Fxq -- '--thread' "$FAKE_LOG"
grep -Fxq -- '3' "$FAKE_LOG"
grep -Fxq -- '--seed' "$FAKE_LOG"
[[ "$(awk -F, 'NR==2 {print $1 ":" $9}' "${WORK}/direct_stats.csv")" == 'astral-pro3:0' ]]

# A STELAR-Pro search-space option must be rejected before any expensive bulk
# simulation starts, with guidance toward the actual ASTRAL-Pro3 controls.
INVALID_DATA_ROOT="${WORK}/must-not-be-created"
if INVALID_LOG=$("${ROOT}/run-bulk-simulated.sh" \
    --method astral-pro3 --simphy-data-dir "$INVALID_DATA_ROOT" \
    --opts-list '--search-space S1' --no-notify 2>&1); then
  echo "Expected the STELAR-Pro-only option to be rejected for ASTRAL-Pro3." >&2
  exit 1
fi
[[ "$INVALID_LOG" == *'--search-space is a STELAR-Pro-only option'* ]]
[[ "$INVALID_LOG" == *'--round N and --subsample N'* ]]
[[ ! -e "$INVALID_DATA_ROOT" ]]

# Drive the public bulk entry point. stat-sim.csv makes sim.sh reuse the fixture.
DATA_ROOT="${WORK}/simphy-data"
RUN_DIR="${DATA_ROOT}/t_4_g_1_sb_0.000001_spmin_100000_spmax_200000/R1"
mkdir -p "$RUN_DIR"
printf '((a,b),(c,d));\n' > "${RUN_DIR}/all_gt.tre"
printf '((a,b),(c,d));\n' > "${RUN_DIR}/s_tree.trees"
printf 'fixture\n' > "${RUN_DIR}/stat-sim.csv"

BULK_LOG=$("${ROOT}/run-bulk-simulated.sh" \
  --method astral-pro3 --taxa-list 4 --genes-list 1 --num-replicates 1 \
  --sb-list 0.000001 --spmin-list 100000 --spmax-list 200000 \
  --simphy-data-dir "$DATA_ROOT" --opts '--thread 3 --seed 7 --verbose 1' \
  --astral-pro3-bin "$FAKE_BIN" --no-gpu-monitor --no-notify 2>&1)
[[ "$BULK_LOG" == *'Method:   astral-pro3'* ]]

RESULTS_DIR="${RUN_DIR}/astral-pro3-outputs/threads_3__seed_7"
[[ -s "${RESULTS_DIR}/out-astral-pro3.tre" ]]
[[ -f "${RESULTS_DIR}/.astral-pro3.success" ]]
[[ -f "${RESULTS_DIR}/.astral-pro3.lock" ]]
[[ "$(awk -F, 'NR==2 {print $1}' "${RESULTS_DIR}/stat-astral-pro3.csv")" == astral-pro3 ]]

SKIP_LOG=$("${ROOT}/test-astral-pro3-simulated.sh" \
  --simphy-data-dir "$DATA_ROOT" -t 4 -g 1 -r R1 \
  --sb 0.000001 --spmin 100000 --spmax 200000 \
  --opts '--thread 3 --seed 7 --verbose 1' \
  --astral-pro3-bin "$FAKE_BIN" --no-gpu-monitor --no-notify 2>&1)
[[ "$SKIP_LOG" == *'SKIPPING: successful output already exists'* ]]

# Notifications use real newlines and the same structured experiment summary as
# STELAR-Pro, including settings, resources, the CSV row, and the stats path.
FAKE_PATH="${WORK}/fake-path"
NOTIFY_LOG="${WORK}/notification.txt"
export NOTIFY_LOG
mkdir -p "$FAKE_PATH"
printf '%s\n' \
  '#!/usr/bin/env bash' \
  'set -euo pipefail' \
  'body=""' \
  'while [[ $# -gt 0 ]]; do' \
  '  case "$1" in -d) body="$2"; shift 2 ;; *) shift ;; esac' \
  'done' \
  'printf "%s" "$body" > "$NOTIFY_LOG"' \
  > "${FAKE_PATH}/curl"
chmod +x "${FAKE_PATH}/curl"
PATH="${FAKE_PATH}:$PATH" NTFY_CHANNEL_NAME=test-channel \
  "${ROOT}/test-astral-pro3-simulated.sh" \
    --simphy-data-dir "$DATA_ROOT" -t 4 -g 1 -r R1 \
    --sb 0.000001 --spmin 100000 --spmax 200000 --fresh \
    --opts '--thread 3 --seed 7 --verbose 1' \
    --astral-pro3-bin "$FAKE_BIN" --no-time-monitor --no-gpu-monitor >/dev/null
grep -q '^✅ ASTRAL-Pro3 completed for 4 taxa and 1 gene trees$' "$NOTIFY_LOG"
grep -q '^Setting: threads_3__seed_7$' "$NOTIFY_LOG"
grep -q '^Options: --thread 3 --seed 7 --verbose 1$' "$NOTIFY_LOG"
grep -q '^RF Rate: 0\.0000$' "$NOTIFY_LOG"
grep -q '^Max CPU RAM: NA MB$' "$NOTIFY_LOG"
grep -q '^alg,setting,num-taxa,gene-trees,replicate,' "$NOTIFY_LOG"
grep -q '^astral-pro3,threads_3__seed_7,4,1,R1,' "$NOTIFY_LOG"
grep -q '^Stats: .*/stat-astral-pro3.csv$' "$NOTIFY_LOG"
if grep -Fq '\n' "$NOTIFY_LOG"; then
  echo 'Notification contains a literal \\n instead of real line breaks.' >&2
  exit 1
fi

"${ROOT}/collect-stats-simulated.sh" \
  --simphy-data-dir "$DATA_ROOT" --out "${WORK}/combined.csv" >/dev/null
grep -q '^astral-pro3,threads_3__seed_7,' "${WORK}/combined.csv"

# The standard-dataset sweep uses the same native wrapper, not the Java ASTRAL
# baseline. "biological" has one configured input, keeping this fixture small.
STANDARD_DATA="${WORK}/standard"
STANDARD_RUN="${STANDARD_DATA}/biological/nuclear/R1"
mkdir -p "$STANDARD_RUN"
printf '((a,b),(c,d));\n' > "${STANDARD_RUN}/all_gt.tre.rooted"
printf '((a,b),(c,d));\n' > "${STANDARD_DATA}/biological/true_tree_trimmed"
"${ROOT}/run-bulk-standard.sh" \
  --base-dir "$WORK" --dataset-dir "$STANDARD_DATA" --stelar-pro-root "$ROOT" \
  --method astral-pro3 --folder biological \
  --opts '--thread 2 --seed 9 --verbose 1' --astral-pro3-bin "$FAKE_BIN" \
  --no-notify >/dev/null
STANDARD_RESULTS="${STANDARD_RUN}/astral-pro3_outputs"
[[ -s "${STANDARD_RESULTS}/output-astral-pro3.tre" ]]
[[ "$(awk -F, 'NR==2 {print $1 ":" $2 ":" $12}' "${STANDARD_RESULTS}/stat-astral-pro3.csv")" == \
    'astral-pro3:threads_2__seed_9:0' ]]
[[ "$(awk -F, 'NR==2 {print $3}' "${STANDARD_RESULTS}/stat-astral-pro3.csv")" == \
    '"--thread 2 --seed 9"' ]]

echo "ASTRAL-Pro3 runner tests: PASS"
