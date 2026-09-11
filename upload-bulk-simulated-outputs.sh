#!/usr/bin/env bash
# Upload the lightweight simulated-results mirror folder by folder.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT}/scripts/phylogeny-data-dir.sh"
source "${ROOT}/scripts/simulated-outputs-mirror.sh"
source "${ROOT}/scripts/hf-python.sh"

OUTPUTS_DIR=""
DATA_DIR=""
SYNC_FIRST=false
METHODS_RAW=""
MIN_TAXA=1
MIN_GENES=1
INCLUDE_INCOMPLETE=true
ALLOW_MISSING_COMMAND=false
REPO_ID="imAniksahA/blab"
REPO_TYPE=dataset
REMOTE_DIR="ph/d/gdl-simulation/outputs"
UPLOADER="${HOME}/utils/hf-data-transfer/hf_upload.py"
PYTHON_BIN=""
DRY_RUN=false
ASSUME_YES=false

usage() {
  cat <<'EOF'
Usage: ./upload-bulk-simulated-outputs.sh [options]

Upload <method>-outputs/<dataset> folders from the results safety copy.

  --outputs-dir DIR            Mirror root (default: $PHYLOGENY_DATA_DIR/outputs/gdl-simulation)
  --simulated-outputs-dir DIR  Alias for --outputs-dir
  --simphy-outputs-dir DIR     Compatibility alias for --outputs-dir
  --sync                       Refresh the mirror before planning
  --data-dir DIR               SimPhy data tree used by --sync
  --methods LIST               Comma/space list (stelar-pro, astral-pro3)
  --min-taxa N                 Minimum taxon count (default: 1)
  --min-gene-trees N           Minimum gene-tree count (default: 1)
  --exclude-incomplete         Skip incomplete datasets
  --allow-missing-command      Permit a dataset without a SimPhy .command
  --repo-id OWNER/REPO         Hugging Face repository
  --repo-type TYPE             dataset, model, or space
  --remote-dir PATH            Remote root (default: ph/d/simulated/outputs)
  --uploader FILE              Folder-aware hf_upload.py
  --python COMMAND             Python with huggingface_hub
  --dry-run                    Validate and print commands only
  --yes, -y                    Skip confirmation
  --help, -h                   Show help
EOF
}

need_value() { [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; exit 2; }; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outputs-dir|--simulated-outputs-dir|--gdl-simulation-outputs-dir|--simphy-outputs-dir) need_value "$@"; OUTPUTS_DIR="$2"; shift 2 ;;
    --sync) SYNC_FIRST=true; shift ;;
    --data-dir|--simphy-data-dir) need_value "$@"; DATA_DIR="$2"; shift 2 ;;
    --methods) need_value "$@"; METHODS_RAW="$2"; shift 2 ;;
    --min-taxa) need_value "$@"; MIN_TAXA="$2"; shift 2 ;;
    --min-gene-trees) need_value "$@"; MIN_GENES="$2"; shift 2 ;;
    --exclude-incomplete) INCLUDE_INCOMPLETE=false; shift ;;
    --allow-missing-command) ALLOW_MISSING_COMMAND=true; shift ;;
    --repo-id) need_value "$@"; REPO_ID="$2"; shift 2 ;;
    --repo-type) need_value "$@"; REPO_TYPE="$2"; shift 2 ;;
    --remote-dir) need_value "$@"; REMOTE_DIR="$2"; shift 2 ;;
    --uploader) need_value "$@"; UPLOADER="$2"; shift 2 ;;
    --python) need_value "$@"; PYTHON_BIN="$2"; shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    --yes|-y) ASSUME_YES=true; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Error: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ "$MIN_TAXA" =~ ^[1-9][0-9]*$ ]] || { echo "Error: --min-taxa must be positive." >&2; exit 2; }
[[ "$MIN_GENES" =~ ^[1-9][0-9]*$ ]] || { echo "Error: --min-gene-trees must be positive." >&2; exit 2; }
[[ "$REPO_ID" =~ ^[^/[:space:]]+/[^/[:space:]]+$ ]] || { echo "Error: invalid --repo-id: $REPO_ID" >&2; exit 2; }
[[ "$REPO_TYPE" =~ ^(dataset|model|space)$ ]] || { echo "Error: invalid --repo-type: $REPO_TYPE" >&2; exit 2; }
REMOTE_DIR="${REMOTE_DIR%/}"
[[ -n "$REMOTE_DIR" && "$REMOTE_DIR" != /* && ! "$REMOTE_DIR" =~ (^|/)[.][.](/|$) ]] || { echo "Error: unsafe --remote-dir: $REMOTE_DIR" >&2; exit 2; }

if [[ "$SYNC_FIRST" == true ]]; then
  DATA_DIR="$(stelar_pro_prepare_simphy_data_dir "$DATA_DIR")" || exit 2
  OUTPUTS_DIR="$(stelar_pro_simulated_outputs_dir "$DATA_DIR" "$OUTPUTS_DIR")" || exit 2
  sync_command=("${ROOT}/sync-simulated-outputs.sh" --simphy-data-dir "$DATA_DIR" --simulated-outputs-dir "$OUTPUTS_DIR" --quiet)
  [[ -n "$METHODS_RAW" ]] && sync_command+=(--methods "$METHODS_RAW")
  [[ "$DRY_RUN" == true ]] && sync_command+=(--dry-run)
  "${sync_command[@]}" || { echo "Error: mirror refresh failed; nothing was uploaded." >&2; exit 1; }
else
  if [[ -z "$OUTPUTS_DIR" ]]; then
    [[ -n "${PHYLOGENY_DATA_DIR:-}" ]] || { echo "Error: PHYLOGENY_DATA_DIR is not set; pass --outputs-dir." >&2; exit 2; }
    OUTPUTS_DIR="${PHYLOGENY_DATA_DIR%/}/outputs/gdl-simulation"
  fi
  OUTPUTS_DIR="$(realpath -m -- "$OUTPUTS_DIR")"
  case "$OUTPUTS_DIR" in /|"${HOME:-/nonexistent}") echo "Error: unsafe outputs root: $OUTPUTS_DIR" >&2; exit 2 ;; esac
  mkdir -p -- "$OUTPUTS_DIR" || exit 2
fi

UPLOADER="$(realpath -m -- "$UPLOADER")"
[[ -f "$UPLOADER" ]] || { echo "Error: uploader not found: $UPLOADER" >&2; exit 2; }
PYTHON_BIN="$(stelar_pro_find_hf_python "$PYTHON_BIN")" || exit 2
if ! "$PYTHON_BIN" "$UPLOADER" --help 2>/dev/null | grep -q -- '--include'; then
  echo "Error: $UPLOADER is not the folder-aware hf_upload.py (its help has no --include)." >&2
  exit 2
fi

declare -A METHOD_FILTER=()
if [[ -n "$METHODS_RAW" ]]; then
  read -r -a method_items <<< "${METHODS_RAW//,/ }"
  for method in "${method_items[@]}"; do
    case "${method,,}" in
      stelar-pro|stelar-pro-outputs|stelar) METHOD_FILTER[stelar-pro-outputs]=1 ;;
      astral-pro3|astral-pro3-outputs|astral-pro|apro3) METHOD_FILTER[astral-pro3-outputs]=1 ;;
      *) echo "Error: unsupported method: $method" >&2; exit 2 ;;
    esac
  done
fi

declare -a SELECTED=()
blocked=0
echo "Simulated outputs uploader"
echo "Outputs directory: $OUTPUTS_DIR"
echo "Destination:       ${REPO_ID}/${REMOTE_DIR}/"
echo
while IFS= read -r -d '' method_path; do
  method_dir="${method_path##*/}"
  [[ ${#METHOD_FILTER[@]} -eq 0 || -n "${METHOD_FILTER[$method_dir]:-}" ]] || continue
  while IFS= read -r -d '' dataset_path; do
    dataset="${dataset_path##*/}"
    stelar_pro_dataset_name_is_valid "$dataset" || continue
    taxa="${BASH_REMATCH[1]}"; genes="${BASH_REMATCH[2]}"; incomplete="${BASH_REMATCH[8]:-}"
    (( taxa >= MIN_TAXA && genes >= MIN_GENES )) || continue
    [[ "$INCLUDE_INCOMPLETE" == true || -z "$incomplete" ]] || continue
    forbidden="$(stelar_pro_find_forbidden_files "$dataset_path" | head -n1)"
    command_file="${dataset_path}/${dataset}.command"
    base_command="${dataset_path}/${dataset%_incomplete}.command"
    result_file="$(find "$dataset_path" -mindepth 3 -type f -print -quit 2>/dev/null)"
    reason=""
    [[ -z "$forbidden" ]] || reason="contains ${forbidden##*/}"
    stelar_pro_mirror_has_unsafe_entries "$dataset_path" && reason="contains symlink/special file"
    [[ -n "$result_file" ]] || reason="contains no result files"
    if [[ "$ALLOW_MISSING_COMMAND" == false && ! -f "$command_file" && ! -f "$base_command" ]]; then reason="missing .command"; fi
    if [[ -n "$reason" ]]; then
      echo "  BLOCKED [$method_dir] $dataset: $reason" >&2
      ((blocked+=1))
    else
      SELECTED+=("${method_dir}/${dataset}")
    fi
  done < <(find "$method_path" -mindepth 1 -maxdepth 1 -type d -name 't_*' -print0 | sort -zV)
done < <(find "$OUTPUTS_DIR" -mindepth 1 -maxdepth 1 -type d -name '*-outputs' -print0 | sort -z)

(( blocked == 0 )) || { echo "Error: blocked unsafe/non-reproducible dataset(s); nothing was uploaded." >&2; exit 1; }
[[ ${#SELECTED[@]} -gt 0 ]] || { echo "No mirrored datasets matched the selection."; exit 0; }

echo "Cases to upload (<method> <dataset> / <replicate> / <setting>):"
for entry in "${SELECTED[@]}"; do
  dataset_path="${OUTPUTS_DIR}/${entry}"
  while IFS= read -r -d '' leaf; do
    relative="${leaf#${dataset_path}/}"
    echo "  [${entry%%/*}] ${entry#*/} / ${relative%%/*} / ${relative#*/}"
  done < <(find "$dataset_path" -mindepth 2 -maxdepth 2 -type d -not -name '.*' -print0 | sort -zV)
done
echo "Plan: ${#SELECTED[@]} dataset folder(s)."

build_upload_command() {
  local entry="$1"
  UPLOAD_COMMAND=("$PYTHON_BIN" "$UPLOADER" --repo-id "$REPO_ID" --repo-type "$REPO_TYPE"
    --local-path "${OUTPUTS_DIR}/${entry}" --path-in-repo "${REMOTE_DIR}/${entry}"
    --commit-message "Simulated outputs: ${entry}")
}

if [[ "$DRY_RUN" == true ]]; then
  for entry in "${SELECTED[@]}"; do build_upload_command "$entry"; printf '  '; printf '%q ' "${UPLOAD_COMMAND[@]}"; printf '\n'; done
  echo "Dry run complete; nothing was uploaded."
  exit 0
fi
if [[ "$ASSUME_YES" == false ]]; then
  read -r -p "Proceed with ${#SELECTED[@]} folder upload(s)? [y/N] " answer
  [[ "$answer" =~ ^[Yy]([Ee][Ss])?$ ]] || { echo "Cancelled; nothing was uploaded."; exit 0; }
fi

failed=0
for entry in "${SELECTED[@]}"; do
  dataset_path="${OUTPUTS_DIR}/${entry}"
  [[ -z "$(stelar_pro_find_forbidden_files "$dataset_path" | head -n1)" ]] || { echo "Error: forbidden file appeared in $entry" >&2; exit 1; }
  stelar_pro_mirror_has_unsafe_entries "$dataset_path" && { echo "Error: unsafe entry appeared in $entry" >&2; exit 1; }
  dataset="${entry#*/}"
  [[ -f "${dataset_path}/${dataset}.command" || -f "${dataset_path}/${dataset%_incomplete}.command" || "$ALLOW_MISSING_COMMAND" == true ]] || {
    echo "Error: command record disappeared from $entry" >&2
    exit 1
  }
  [[ -n "$(find "$dataset_path" -mindepth 3 -type f -print -quit 2>/dev/null)" ]] || {
    echo "Error: result files disappeared from $entry" >&2
    exit 1
  }
done
for entry in "${SELECTED[@]}"; do
  build_upload_command "$entry"
  echo "Uploading $entry"
  "${UPLOAD_COMMAND[@]}" || ((failed+=1))
done
(( failed == 0 )) || exit 1
echo "Uploaded ${#SELECTED[@]} dataset folder(s)."
