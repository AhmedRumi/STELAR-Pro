#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "${ROOT}/scripts/simulated-outputs-mirror.sh"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelar-pro-mirror-test.XXXXXX")"
trap 'status=$?; rm -rf -- "$WORK"; exit "$status"' EXIT

DATA="${WORK}/simphy/data"
OUTPUTS="${WORK}/outputs/gdl-simulation"
DATASET="t_10_g_20_sb_0.000001_spmin_100000_spmax_200000"
STELAR_SOURCE="${DATA}/${DATASET}/R1/stelar-pro-outputs/search-space_S1"
ASTRAL_SOURCE="${DATA}/${DATASET}/R1/astral-pro3-outputs/default"
mkdir -p "$STELAR_SOURCE" "$ASTRAL_SOURCE"

[[ "$(stelar_pro_simulated_outputs_dir "$DATA" "")" == "$OUTPUTS" ]]
if stelar_pro_simulated_outputs_dir "$DATA" "${DATA}/mirror" >/dev/null 2>&1; then
  echo "Expected an in-data mirror root to be rejected." >&2
  exit 1
fi

printf 'sim command\n' > "${DATA}/${DATASET}/${DATASET}.command"
printf 'sim params\n' > "${DATA}/${DATASET}/${DATASET}.params"
printf 'database\n' > "${DATA}/${DATASET}/${DATASET}.db"
printf 'gene trees\n' > "${DATA}/${DATASET}/R1/all_gt.tre"
printf 'result\n' > "${STELAR_SOURCE}/out-stelar-pro.tre"
printf 'stats\n' > "${STELAR_SOURCE}/stat-stelar-pro.csv"
printf 'command\n' > "${STELAR_SOURCE}/out-stelar-pro.command"
touch "${STELAR_SOURCE}/.stelar-pro.success"

stelar_pro_mirror_simulated_results "$DATA" "$OUTPUTS" "$STELAR_SOURCE" >/dev/null
STELAR_DEST="${OUTPUTS}/stelar-pro-outputs/${DATASET}/R1/search-space_S1"
[[ -s "${STELAR_DEST}/out-stelar-pro.tre" ]]
[[ -f "${STELAR_DEST}/.stelar-pro.success" ]]
[[ -s "${OUTPUTS}/stelar-pro-outputs/${DATASET}/${DATASET}.command" ]]
[[ -s "${OUTPUTS}/stelar-pro-outputs/${DATASET}/${DATASET}.params" ]]
[[ ! -e "${OUTPUTS}/stelar-pro-outputs/${DATASET}/${DATASET}.db" ]]
[[ ! -e "${STELAR_DEST}/all_gt.tre" ]]

# A refresh replaces the exact result leaf, so files removed at the source do
# not linger in the safety copy.
printf 'stale\n' > "${STELAR_DEST}/stale.txt"
stelar_pro_mirror_simulated_results "$DATA" "$OUTPUTS" "$STELAR_SOURCE" >/dev/null
[[ ! -e "${STELAR_DEST}/stale.txt" ]]

printf 'astral result\n' > "${ASTRAL_SOURCE}/out-astral-pro3.tre"
touch "${ASTRAL_SOURCE}/.astral-pro3.success"
stelar_pro_mirror_simulated_results "$DATA" "$OUTPUTS" "$ASTRAL_SOURCE" >/dev/null
[[ -s "${OUTPUTS}/astral-pro3-outputs/${DATASET}/R1/default/out-astral-pro3.tre" ]]
[[ -f "${OUTPUTS}/astral-pro3-outputs/${DATASET}/R1/default/.astral-pro3.success" ]]

# A forbidden source file aborts before replacing a previously safe copy.
printf 'do not copy\n' > "${ASTRAL_SOURCE}/all_gt.tre"
if stelar_pro_mirror_simulated_results "$DATA" "$OUTPUTS" "$ASTRAL_SOURCE" >/dev/null 2>&1; then
  echo "Expected a result leaf containing gene trees to be rejected." >&2
  exit 1
fi
[[ -s "${OUTPUTS}/astral-pro3-outputs/${DATASET}/R1/default/out-astral-pro3.tre" ]]
[[ ! -e "${OUTPUTS}/astral-pro3-outputs/${DATASET}/R1/default/all_gt.tre" ]]
rm -f "${ASTRAL_SOURCE}/all_gt.tre"

rm -rf "$OUTPUTS"
"${ROOT}/sync-simulated-outputs.sh" --simphy-data-dir "$DATA" --simulated-outputs-dir "$OUTPUTS" --quiet >/dev/null
[[ -s "${STELAR_DEST}/out-stelar-pro.tre" ]]
[[ -s "${OUTPUTS}/astral-pro3-outputs/${DATASET}/R1/default/out-astral-pro3.tre" ]]

# Exercise uploader discovery without network access. The fake Python validates
# interpreter detection and exposes the folder-upload flag in helper --help.
FAKE_PYTHON="${WORK}/python-with-hf"
FAKE_UPLOADER="${WORK}/hf_upload.py"
printf '# fixture\n' > "$FAKE_UPLOADER"
printf '%s\n' '#!/usr/bin/env bash' 'if [[ "${1:-}" == -c ]]; then exit 0; fi' 'echo --include' > "$FAKE_PYTHON"
chmod +x "$FAKE_PYTHON"
UPLOAD_PLAN="$("${ROOT}/upload-bulk-simulated-outputs.sh" --outputs-dir "$OUTPUTS" \
  --methods 'stelar-pro,astral-pro3' --python "$FAKE_PYTHON" --uploader "$FAKE_UPLOADER" --dry-run)"
[[ "$UPLOAD_PLAN" == *"[stelar-pro-outputs] ${DATASET} / R1 / search-space_S1"* ]]
[[ "$UPLOAD_PLAN" == *"[astral-pro3-outputs] ${DATASET} / R1 / default"* ]]
[[ "$UPLOAD_PLAN" == *'Dry run complete; nothing was uploaded.'* ]]

echo "Simulated outputs mirror tests: PASS"
