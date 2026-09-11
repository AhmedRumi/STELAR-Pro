#!/usr/bin/env bash
# Shared, side-effect-free helpers for the lightweight simulated-results mirror.

stelar_pro_print_shell_command() {
  local argument first=true
  for argument in "$@"; do
    if [[ "$first" == true ]]; then first=false; else printf ' '; fi
    printf '%q' "$argument"
  done
  printf '\n'
}

stelar_pro_simulated_outputs_dir() {
  local data_dir="$1" requested="${2:-}" resolved base
  data_dir="$(realpath -m -- "$data_dir")"

  if [[ -n "$requested" ]]; then
    if [[ "$requested" == "~/"* ]]; then
      [[ -n "${HOME:-}" ]] || { echo "Error: cannot expand '~' because HOME is not set." >&2; return 2; }
      requested="${HOME}/${requested:2}"
    fi
    resolved="$(realpath -m -- "$requested")"
  elif [[ "$data_dir" == */simphy/data ]]; then
    base="${data_dir%/simphy/data}"
    resolved="${base}/outputs/gdl-simulation"
  elif [[ "$(basename "$data_dir")" == data ]]; then
    resolved="$(dirname "$data_dir")/outputs/gdl-simulation"
  else
    resolved="${data_dir}_outputs/gdl-simulation"
  fi

  case "${resolved}/" in
    "${data_dir}/"*)
      echo "Error: outputs mirror cannot be inside the SimPhy data directory: $resolved" >&2
      return 2
      ;;
  esac
  case "${data_dir}/" in
    "${resolved}/"*)
      echo "Error: outputs mirror cannot contain the SimPhy data directory: $resolved" >&2
      return 2
      ;;
  esac
  case "$resolved" in
    /|"${HOME:-/nonexistent}")
      echo "Error: refusing unsafe outputs mirror root: $resolved" >&2
      return 2
      ;;
  esac

  printf '%s\n' "$resolved"
}

stelar_pro_dataset_name_is_valid() {
  [[ "$1" =~ ^t_([1-9][0-9]*)_g_([1-9][0-9]*)_sb_([0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?)_spmin_([1-9][0-9]*)_spmax_([1-9][0-9]*)(_incomplete)?$ ]]
}

stelar_pro_mirror_has_forbidden_files() {
  local directory="$1"
  [[ -n "$(stelar_pro_find_forbidden_files "$directory" | head -n1)" ]]
}

stelar_pro_find_forbidden_files() {
  local directory="$1"
  find "$directory" -type f \
    \( -name all_gt.tre -o -name s_tree.trees -o -name l_trees.trees \
       -o -name 'g_trees*.trees' -o -name '*.db' -o -name '*.db-journal' \
       -o -name '*.zip' -o -name stat-sim.csv \) -print 2>/dev/null
}

stelar_pro_mirror_has_unsafe_entries() {
  local directory="$1"
  [[ -n "$(find "$directory" \( -type l -o \( ! -type d ! -type f \) \) -print -quit 2>/dev/null)" ]]
}

stelar_pro_copy_dataset_records() {
  local data_dir="$1" outputs_root="$2" dataset="$3" method_dir="$4"
  local source_dataset="${data_dir}/${dataset}"
  local destination="${outputs_root}/${method_dir}/${dataset}"
  local base_dataset="$dataset" record target temporary copied_command=false
  [[ "$dataset" == *_incomplete ]] && base_dataset="${dataset%_incomplete}"

  mkdir -p -- "$destination"
  for source_dataset in "${data_dir}/${base_dataset}" "${data_dir}/${dataset}"; do
    [[ -d "$source_dataset" ]] || continue
    while IFS= read -r -d '' record; do
      [[ "$record" == *.command ]] && copied_command=true
      target="${destination}/$(basename "$record")"
      if [[ -f "$target" ]] && cmp -s -- "$record" "$target"; then
        continue
      fi
      temporary="${target}.tmp.$$"
      if cp -p -- "$record" "$temporary" && mv -f -- "$temporary" "$target"; then
        :
      else
        rm -f -- "$temporary" 2>/dev/null || true
        return 1
      fi
    done < <(find "$source_dataset" -maxdepth 1 -type f \
      \( -name '*.command' -o -name '*.params' \) -print0 | sort -z)
  done

  if [[ "$copied_command" == false ]]; then
    echo "Warning: no SimPhy .command record found for dataset '$dataset'." >&2
  fi
}

# Mirror a single result leaf:
#   DATA/<dataset>/<R>/<method>-outputs/<setting>
# becomes:
#   OUTPUTS/<method>-outputs/<dataset>/<R>/<setting>
stelar_pro_mirror_simulated_results() {
  local data_dir outputs_root results_dir relative dataset replicate method_dir setting extra
  local destination_parent destination temporary backup=""
  data_dir="$(realpath -m -- "$1")"
  outputs_root="$(stelar_pro_simulated_outputs_dir "$data_dir" "$2")" || return
  results_dir="$(realpath -m -- "$3")"

  [[ -d "$results_dir" ]] || {
    echo "Error: results directory does not exist: $results_dir" >&2
    return 2
  }
  case "${results_dir}/" in
    "${data_dir}/"*) ;;
    *) echo "Error: results directory is outside the SimPhy data tree: $results_dir" >&2; return 2 ;;
  esac

  relative="${results_dir#${data_dir}/}"
  IFS=/ read -r dataset replicate method_dir setting extra <<< "$relative"
  if [[ -n "${extra:-}" || -z "$setting" ]] || \
      ! stelar_pro_dataset_name_is_valid "$dataset" || \
      [[ ! "$replicate" =~ ^R[1-9][0-9]*$ ]] || \
      [[ ! "$method_dir" =~ ^[a-z0-9][a-z0-9-]*-outputs$ ]] || \
      [[ "$setting" == . || "$setting" == .. ]]; then
    echo "Error: results directory does not have the required dataset/R/method-outputs/setting shape: $results_dir" >&2
    return 2
  fi
  if stelar_pro_mirror_has_forbidden_files "$results_dir"; then
    echo "Error: refusing to mirror simulated input/archive files from: $results_dir" >&2
    return 2
  fi
  if stelar_pro_mirror_has_unsafe_entries "$results_dir"; then
    echo "Error: refusing to mirror symlinks or special files from: $results_dir" >&2
    return 2
  fi

  stelar_pro_copy_dataset_records "$data_dir" "$outputs_root" "$dataset" "$method_dir" || return 1

  destination_parent="${outputs_root}/${method_dir}/${dataset}/${replicate}"
  destination="${destination_parent}/${setting}"
  mkdir -p -- "$destination_parent"
  temporary="$(mktemp -d "${destination_parent}/.${setting}.mirror.XXXXXX")" || return 1
  if ! cp -a -- "${results_dir}/." "$temporary/"; then
    rm -rf -- "$temporary"
    return 1
  fi
  if stelar_pro_mirror_has_forbidden_files "$temporary"; then
    echo "Error: forbidden file appeared while mirroring: $results_dir" >&2
    rm -rf -- "$temporary"
    return 2
  fi

  if [[ -e "$destination" ]]; then
    backup="$(mktemp -d "${destination_parent}/.${setting}.previous.XXXXXX")" || {
      rm -rf -- "$temporary"
      return 1
    }
    rmdir -- "$backup" || {
      rm -rf -- "$temporary" "$backup"
      return 1
    }
    if ! mv -- "$destination" "$backup"; then
      rm -rf -- "$temporary"
      return 1
    fi
  fi
  if ! mv -- "$temporary" "$destination"; then
    [[ -n "$backup" && -e "$backup" ]] && mv -- "$backup" "$destination" 2>/dev/null || true
    rm -rf -- "$temporary"
    return 1
  fi
  [[ -n "$backup" && -e "$backup" ]] && rm -rf -- "$backup"

  printf '%s\n' "$destination"
}

stelar_pro_git_revision() {
  local root="$1" revision=unknown
  if git -C "$root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    revision="$(git -C "$root" rev-parse HEAD 2>/dev/null || echo unknown)"
    [[ -n "$(git -C "$root" status --porcelain 2>/dev/null)" ]] && revision+="-dirty"
  fi
  printf '%s\n' "$revision"
}
