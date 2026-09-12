#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelar-gdl-test.XXXXXX")"
trap 'rm -rf -- "$WORK"' EXIT
DATA="${WORK}/gdl-simulation/data"
OUTPUTS="${WORK}/outputs/gdl-simulation"
DATASET=taxa3_gt2_dup1_loss0_pop50000000
RUN="${DATA}/${DATASET}/R1"
mkdir -p "${RUN}/species-tree" "${RUN}/true-genetrees" "${DATA}/${DATASET}/simphy_raw"
printf '((0,1),(2,3));\n' >"${RUN}/species-tree/s_tree.trees"
printf '(((0_0_0,0_1_0),1_0_0),(2_0_0,3_0_0));\n((0_0_0,1_0_0),(2_0_0,3_0_0));\n' >"${RUN}/true-genetrees/all_gt.trees"
printf 'parameters\n' >"${DATA}/${DATASET}/params.txt"
printf 'generator command\n' >"${DATA}/${DATASET}/simphy_raw/simphy_raw.command"
printf 'generator params\n' >"${DATA}/${DATASET}/simphy_raw/simphy_raw.params"
cp "${ROOT}/test/fixtures/fake_gdl_astral.py" "${WORK}/astral-pro3"
chmod +x "${WORK}/astral-pro3"
export GDL_CALLS="${WORK}/calls.log"
COMMON=(--gdl-data-dir "$DATA" --taxa-list 3 --gt-list 2 --method astral-pro3 --opts '--thread 2'
  --astral-pro3-bin "${WORK}/astral-pro3" --no-notify --yes --no-time-monitor)
PLAN="$("${ROOT}/run-bulk-gdl-simulated.sh" "${COMMON[@]}" --dry-run)"
[[ "$PLAN" == *'true-genetrees/all_gt.trees'* && "$PLAN" == *'species-tree/s_tree.trees'* ]]
[[ ! -e "$GDL_CALLS" && ! -e "$OUTPUTS" && ! -e "${RUN}/astral-pro3-outputs" ]]
# Missing strict selection must abort before any method starts.
if "${ROOT}/run-bulk-gdl-simulated.sh" "${COMMON[@]}" --num-replicates 2 --strict-missing >"${WORK}/strict.log" 2>&1; then exit 1; fi
[[ ! -e "$GDL_CALLS" ]]
"${ROOT}/run-bulk-gdl-simulated.sh" "${COMMON[@]}" --num-replicates 2 >"${WORK}/bulk.log" 2>&1
RESULT="${RUN}/astral-pro3-outputs/threads_2"
DEST="${OUTPUTS}/astral-pro3-outputs/${DATASET}/R1/threads_2"
[[ -s "${RESULT}/out-astral-pro3.tre" && -f "${RESULT}/.astral-pro3.success" && -f "${RESULT}/.astral-pro3.lock" ]]
[[ -s "${DEST}/input-preparation.json" && -s "${DEST}/out-astral-pro3.command" ]]
[[ -s "${OUTPUTS}/astral-pro3-outputs/${DATASET}/simphy_raw.command" && -s "${OUTPUTS}/astral-pro3-outputs/${DATASET}/params.txt" ]]
[[ "$(awk -F, 'NR==2 {print $9":"$13":"$15":"$16}' "${RESULT}/stat-astral-pro3.csv")" == '0.0000000000:NA:0:ok' ]]
[[ "$(wc -l <"$GDL_CALLS")" == 1 ]]
"${ROOT}/run-bulk-gdl-simulated.sh" "${COMMON[@]}" >"${WORK}/skip.log" 2>&1
[[ "$(wc -l <"$GDL_CALLS")" == 1 ]]
rg -q 'SKIPPING' "${WORK}/skip.log"
# A failed fresh run must not overwrite either successful original or mirror.
if GDL_FAIL=1 "${ROOT}/run-bulk-gdl-simulated.sh" "${COMMON[@]}" --fresh >"${WORK}/fail.log" 2>&1; then exit 1; fi
[[ -f "${RESULT}/.astral-pro3.success" && -f "${DEST}/.astral-pro3.success" ]]
mapfile -t FAILED < <(find "${RUN}/astral-pro3-outputs" -maxdepth 1 -type d -name 'threads_2__failed_*')
[[ ${#FAILED[@]} == 1 && ! -e "${FAILED[0]}/.astral-pro3.success" && ! -e "${FAILED[0]}/.astral-pro3.lock" ]]
"${ROOT}/run-bulk-gdl-simulated.sh" "${COMMON[@]}" --fresh >"${WORK}/fresh.log" 2>&1
[[ -n "$(find "${RUN}/astral-pro3-outputs" -maxdepth 1 -type d -name 'threads_2__previous_*')" ]]
# Fingerprint change must trigger a safe rerun rather than stale skip.
printf '\n' >>"${RUN}/true-genetrees/all_gt.trees"
COUNT="$(wc -l <"$GDL_CALLS")"
"${ROOT}/run-bulk-gdl-simulated.sh" "${COMMON[@]}" >"${WORK}/changed.log" 2>&1
[[ "$(wc -l <"$GDL_CALLS")" == "$((COUNT + 1))" ]]
# Collector includes only current successful settings, never rerun archives.
"${ROOT}/collect-stats-simulated.sh" --gdl-data-dir "$DATA" --out "${WORK}/combined.csv" >"${WORK}/collect.log"
[[ "$(wc -l <"${WORK}/combined.csv")" == 2 ]]
rg -q 'replicate,dup,loss,pop' "${WORK}/combined.csv"
"${ROOT}/sync-simulated-outputs.sh" --gdl-data-dir "$DATA" --quiet >"${WORK}/sync.log"
source "${ROOT}/scripts/simulated-outputs-mirror.sh"
[[ -z "$(stelar_pro_find_forbidden_files "$OUTPUTS")" ]]
[[ -z "$(find "$OUTPUTS" -name species-labelled.trees -o -name separate)" ]]
# Existing uploader must recognize GDL generator provenance without network.
printf '#!/usr/bin/env bash\nif [[ "${1:-}" == -c ]]; then exit 0; fi\necho --include\n' >"${WORK}/python-hf"
chmod +x "${WORK}/python-hf"
printf '# fixture\n' >"${WORK}/uploader.py"
"${ROOT}/upload-bulk-simulated-outputs.sh" --outputs-dir "$OUTPUTS" --python "${WORK}/python-hf" --uploader "${WORK}/uploader.py" --dry-run >"${WORK}/upload.log"
rg -q "${DATASET} / R1 / threads_2" "${WORK}/upload.log"
# Complete label normalization and counts are required, even when native opts fail.
if "${ROOT}/run-bulk-gdl-simulated.sh" "${COMMON[@]}" --opts '--search-space S1' >"${WORK}/invalid.log" 2>&1; then exit 1; fi
if GDL_EMPTY=1 "${ROOT}/run-bulk-gdl-simulated.sh" "${COMMON[@]}" --fresh >"${WORK}/empty.log" 2>&1; then exit 1; fi
# Capture actual notification data locally, without accessing ntfy.
mkdir -p "${WORK}/fake-path"
export GDL_NOTIFICATION="${WORK}/notification.txt"
printf '%s\n' '#!/usr/bin/env bash' 'while (($#)); do' 'if [[ "$1" == -d ]]; then printf "%s" "$2" >"$GDL_NOTIFICATION"; exit 0; fi' 'shift' 'done' >"${WORK}/fake-path/curl"
chmod +x "${WORK}/fake-path/curl"
PATH="${WORK}/fake-path:${PATH}" "${ROOT}/test-gdl-simulated.sh" --gdl-data-dir "$DATA" --dataset "$DATASET" \
  --method astral-pro3 --opts '--thread 3' --astral-pro3-bin "${WORK}/astral-pro3" --no-time-monitor >"${WORK}/notify.log" 2>&1
[[ -s "$GDL_NOTIFICATION" && "$(wc -l <"$GDL_NOTIFICATION")" -gt 10 ]]
rg -q '^Dataset: taxa3_gt2_dup1_loss0_pop50000000$' "$GDL_NOTIFICATION"
rg -q '^Replicate: R1$' "$GDL_NOTIFICATION"
rg -q '^Options: --thread 3$' "$GDL_NOTIFICATION"
rg -q '^RF Rate: 0.0000000000 \(ok\)$' "$GDL_NOTIFICATION"
rg -q '^Max GPU VRAM: NA MB$' "$GDL_NOTIFICATION"
rg -q 'replicate,dup,loss,pop' "$GDL_NOTIFICATION"
# GDL combined input is explicitly forbidden even if placed in a result leaf.
printf 'raw genes\n' >"${RESULT}/all_gt.trees"
if stelar_pro_mirror_simulated_results "$DATA" "$OUTPUTS" "$RESULT" >"${WORK}/forbidden.log" 2>&1; then exit 1; fi
[[ ! -e "${DEST}/all_gt.trees" ]]
echo 'GDL simulated runner tests: PASS'
