#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/stelar-pro-success-test.XXXXXX")"
MIRROR_ROOT="${TMP}-outputs/gdl-simulation"
trap 'status=$?; rm -rf -- "$TMP" "${TMP}-outputs"; exit "$status"' EXIT

DATASET_NAME="t_4_g_1_sb_0.000001_spmin_100000_spmax_200000"
RUN_DIR="$TMP/$DATASET_NAME/R1"
mkdir -p "$RUN_DIR"
printf 'simphy fixture command\n' > "$TMP/$DATASET_NAME/$DATASET_NAME.command"
printf 'simphy fixture params\n' > "$TMP/$DATASET_NAME/$DATASET_NAME.params"
printf '((a,b),(c,d));\n' > "$RUN_DIR/all_gt.tre"
printf '((a,b),(c,d));\n' > "$RUN_DIR/s_tree.trees"

COMMON=(--simphy-data-dir "$TMP" -t 4 -g 1 -r R1
  --simulated-outputs-dir "$MIRROR_ROOT"
  --sb 0.000001 --spmin 100000 --spmax 200000
  --opts '--search-space S1 --cpu -q'
  --no-time-monitor --no-gpu-monitor --no-notify)

"$ROOT/test-stelar-pro-simulated.sh" "${COMMON[@]}" >/dev/null
RESULTS_DIR=$(find "$RUN_DIR/stelar-pro-outputs" -mindepth 1 -maxdepth 1 -type d)
OUTPUT="$RESULTS_DIR/out-stelar-pro.tre"
SIDE="$RESULTS_DIR/out-stelar-pro_stats.csv"
SUCCESS="$RESULTS_DIR/.stelar-pro.success"
[[ -s "$OUTPUT" && -s "$SIDE" && -s "$SUCCESS" ]]
MIRRORED_RESULTS="$MIRROR_ROOT/stelar-pro-outputs/$DATASET_NAME/R1/$(basename "$RESULTS_DIR")"
[[ -s "$MIRRORED_RESULTS/out-stelar-pro.tre" ]]
[[ -s "$MIRRORED_RESULTS/out-stelar-pro.command" ]]
[[ -f "$MIRRORED_RESULTS/.stelar-pro.success" ]]
[[ -s "$MIRROR_ROOT/stelar-pro-outputs/$DATASET_NAME/$DATASET_NAME.command" ]]
STAT_FILE="$RESULTS_DIR/stat-stelar-pro.csv"
grep -q 'optimal-triplet-score' "$STAT_FILE"
! grep -qi 'quartet' "$STAT_FILE"
[[ "$(awk -F, 'NR==2 {print $10}' "$STAT_FILE")" =~ ^[0-9]+([.][0-9]+)?$ ]]

# A failed sidecar plus a stale tree must never be accepted as completed.
rm -f "$SUCCESS"
sed -i '2s/,0$/,1/' "$SIDE"
rerun_log=$("$ROOT/test-stelar-pro-simulated.sh" "${COMMON[@]}" 2>&1)
[[ "$rerun_log" == *"Previous statistics exist but no successful output was recorded; rerunning."* ]]
[[ -s "$OUTPUT" && -s "$SUCCESS" ]]
[[ "$(awk -F, 'NR==2 {print $9}' "$SIDE")" == "0" ]]

printf 'stale\n' > "$MIRRORED_RESULTS/stale-after-run.txt"
skip_log=$("$ROOT/test-stelar-pro-simulated.sh" "${COMMON[@]}" 2>&1)
[[ "$skip_log" == *"SKIPPING: successful output already exists"* ]]
[[ ! -e "$MIRRORED_RESULTS/stale-after-run.txt" ]]
diff -qr "$RESULTS_DIR" "$MIRRORED_RESULTS" >/dev/null

COMBINED="${TMP}/combined.csv"
"${ROOT}/collect-stats-simulated.sh" --simphy-data-dir "$TMP" --out "$COMBINED" >/dev/null
grep -q 'optimal-triplet-score' "$COMBINED"
! grep -qi 'quartet' "$COMBINED"
[[ "$(awk -F, 'NR==2 {print $10}' "$COMBINED")" =~ ^[0-9]+([.][0-9]+)?$ ]]

echo "Simulated-run success detection: PASS"
