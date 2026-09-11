#!/usr/bin/env bash
# Backfill the lightweight, method-first mirror of existing simulated results.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT}/scripts/phylogeny-data-dir.sh"
source "${ROOT}/scripts/simulated-outputs-mirror.sh"

SIMPHY_DATA_DIR=""
SIMULATED_OUTPUTS_DIR=""
METHODS_RAW="stelar-pro,astral-pro3"
DRY_RUN=false
QUIET=false

usage() {
  cat <<'EOF'
Usage: ./sync-simulated-outputs.sh [options]

Copy existing method result leaves into the canonical method-first mirror.
Simulation inputs, databases, archives, and stat-sim.csv are never copied.

Options:
  --simphy-data-dir DIR         SimPhy data tree
  --simulated-outputs-dir DIR   Mirror root
                                (default: $PHYLOGENY_DATA_DIR/outputs/gdl-simulation)
  --gdl-simulation-outputs-dir DIR
                                Alias for --simulated-outputs-dir
  --simphy-outputs-dir DIR      Compatibility alias
  --methods LIST                Comma/space list: stelar-pro, astral-pro3
  --dry-run                     Show leaves without copying
  --quiet                       Print only the final count and errors
  --help, -h                    Show help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --simphy-data-dir|--data-dir) SIMPHY_DATA_DIR="$2"; shift 2 ;;
    --simulated-outputs-dir|--gdl-simulation-outputs-dir|--simphy-outputs-dir) SIMULATED_OUTPUTS_DIR="$2"; shift 2 ;;
    --methods) METHODS_RAW="$2"; shift 2 ;;
    --dry-run) DRY_RUN=true; shift ;;
    --quiet) QUIET=true; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Error: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

SIMPHY_DATA_DIR="$(stelar_pro_prepare_simphy_data_dir "$SIMPHY_DATA_DIR")"
SIMULATED_OUTPUTS_DIR="$(stelar_pro_simulated_outputs_dir "$SIMPHY_DATA_DIR" "$SIMULATED_OUTPUTS_DIR")"

declare -A WANTED_METHODS=()
read -r -a requested_methods <<< "${METHODS_RAW//,/ }"
for method in "${requested_methods[@]}"; do
  case "${method,,}" in
    stelar-pro|stelar-pro-outputs|stelar) WANTED_METHODS[stelar-pro-outputs]=1 ;;
    astral-pro3|astral-pro3-outputs|astral-pro|apro3) WANTED_METHODS[astral-pro3-outputs]=1 ;;
    *) echo "Error: unsupported method in --methods: $method" >&2; exit 2 ;;
  esac
done
[[ ${#WANTED_METHODS[@]} -gt 0 ]] || { echo "Error: --methods cannot be empty." >&2; exit 2; }

count=0
failed=0
while IFS= read -r -d '' results_dir; do
  relative="${results_dir#${SIMPHY_DATA_DIR}/}"
  IFS=/ read -r dataset replicate method_dir setting extra <<< "$relative"
  [[ -z "${extra:-}" && -n "${setting:-}" ]] || continue
  [[ -n "${WANTED_METHODS[$method_dir]:-}" ]] || continue
  stelar_pro_dataset_name_is_valid "$dataset" || continue
  [[ "$replicate" =~ ^R[1-9][0-9]*$ ]] || continue

  if stelar_pro_mirror_has_forbidden_files "$results_dir"; then
    echo "Error: refusing result leaf containing simulated input/archive files: $results_dir" >&2
    ((failed+=1))
    continue
  fi
  if stelar_pro_mirror_has_unsafe_entries "$results_dir"; then
    echo "Error: refusing result leaf containing symlinks or special files: $results_dir" >&2
    ((failed+=1))
    continue
  fi

  destination="${SIMULATED_OUTPUTS_DIR}/${method_dir}/${dataset}/${replicate}/${setting}"
  if [[ "$QUIET" == false ]]; then
    printf '%s%s -> %s\n' "$(if [[ "$DRY_RUN" == true ]]; then printf '[dry-run] '; fi)" "$results_dir" "$destination"
  fi
  if [[ "$DRY_RUN" == false ]]; then
    stelar_pro_mirror_simulated_results "$SIMPHY_DATA_DIR" "$SIMULATED_OUTPUTS_DIR" "$results_dir" >/dev/null
  fi
  ((count+=1)) || true
done < <(find "$SIMPHY_DATA_DIR" -mindepth 4 -maxdepth 4 -type d -print0 | sort -z)

if (( failed > 0 )); then
  echo "Failed safety checks for $failed result leaf/leaves; the safe leaves were processed." >&2
  exit 1
fi

if [[ "$DRY_RUN" == true ]]; then
  echo "Would mirror $count result leaf/leaves into $SIMULATED_OUTPUTS_DIR."
else
  echo "Mirrored $count result leaf/leaves into $SIMULATED_OUTPUTS_DIR."
fi
