#!/usr/bin/env bash

stelar_pro_find_hf_python() {
  local requested="${1:-}" candidate resolved tried=""
  local -a candidates=()

  if [[ -n "$requested" ]]; then
    if [[ "$requested" == */* ]]; then resolved="$requested"; else resolved="$(command -v "$requested" 2>/dev/null || true)"; fi
    [[ -n "$resolved" && -x "$resolved" ]] || { echo "Error: Python interpreter is not executable: $requested" >&2; return 2; }
    "$resolved" -c 'import huggingface_hub' >/dev/null 2>&1 || {
      echo "Error: $resolved cannot import huggingface_hub." >&2
      return 2
    }
    printf '%s\n' "$resolved"
    return
  fi

  candidates=(python3 python)
  [[ -n "${CONDA_EXE:-}" ]] && candidates+=("${CONDA_EXE%/bin/conda}/bin/python")
  [[ -n "${CONDA_PREFIX:-}" ]] && candidates+=("${CONDA_PREFIX}/bin/python")
  candidates+=("${HOME}/miniconda3/bin/python" "${HOME}/anaconda3/bin/python" /usr/bin/python3)
  for candidate in "${candidates[@]}"; do
    if [[ "$candidate" == */* ]]; then resolved="$candidate"; else resolved="$(command -v "$candidate" 2>/dev/null || true)"; fi
    [[ -n "$resolved" && -x "$resolved" ]] || continue
    tried+="${tried:+, }$resolved"
    if "$resolved" -c 'import huggingface_hub' >/dev/null 2>&1; then printf '%s\n' "$resolved"; return; fi
  done
  echo "Error: no Python interpreter with huggingface_hub was found (tried: ${tried:-none})." >&2
  echo "Install huggingface_hub/hf_xet or pass --python /path/to/python." >&2
  return 2
}
