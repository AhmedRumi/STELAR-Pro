#!/usr/bin/env bash
# Shared setting-name encoder for the experiment runners.

sanitize_setting_part() {
  local value="$1"
  value="${value// /-}"
  value="${value//\//-}"
  value="${value//:/-}"
  value="${value//=/-}"
  value="${value//,/.-}"
  printf '%s' "$value"
}

canonical_search_space_name() {
  local value="${1,,}"
  value="${value//_/-}"
  case "$value" in
    1|s1) printf 'S1' ;;
    2|s2) printf 'S2' ;;
    3|s3) printf 'S3' ;;
    *) sanitize_setting_part "$1" ;;
  esac
}

# Encode meaningful options as option_value groups separated by '__'.
# Display-only verbosity flags do not describe an experimental setting and are
# omitted. Friendly and legacy STELAR-Pro aliases use the compact preset names.
build_setting_name_from_opts() {
  local raw="$1"
  local -a tokens=()
  local -a parts=()
  local i=0 token key value next

  if [[ -z "${raw// }" ]]; then
    printf 'default'
    return
  fi

  read -r -a tokens <<< "$raw"
  while (( i < ${#tokens[@]} )); do
    token="${tokens[$i]}"
    case "$token" in
      -v|-vv|-vvv|-q|--quiet|--verbose)
        ((i+=1))
        continue
        ;;
      --intersection-method|--im|--weight-intersection-method)
        echo "Error: $token was removed from STELAR-Pro." >&2
        return 2
        ;;
      --search-space)
        if (( i + 1 >= ${#tokens[@]} )); then
          ((i+=1))
          continue
        fi
        value="${tokens[$((i + 1))]}"
        parts+=("search-space_$(canonical_search_space_name "$value")")
        ((i+=2))
        ;;
      --intersection-method=*|--im=*|--weight-intersection-method=*)
        echo "Error: ${token%%=*} was removed from STELAR-Pro." >&2
        return 2
        ;;
      --search-space=*)
        value="${token#*=}"
        parts+=("search-space_$(canonical_search_space_name "$value")")
        ((i+=1))
        ;;
      -t|-T|--threads|--num-threads)
        if (( i + 1 < ${#tokens[@]} )); then
          parts+=("threads_$(sanitize_setting_part "${tokens[$((i + 1))]}")")
          ((i+=2))
        else
          ((i+=1))
        fi
        ;;
      -m|--seeds)
        if (( i + 1 < ${#tokens[@]} )); then
          parts+=("seeds_$(sanitize_setting_part "${tokens[$((i + 1))]}")")
          ((i+=2))
        else
          ((i+=1))
        fi
        ;;
      --*=*)
        key="${token%%=*}"
        value="${token#*=}"
        parts+=("$(sanitize_setting_part "${key#--}")_$(sanitize_setting_part "$value")")
        ((i+=1))
        ;;
      --*)
        key="${token#--}"
        if (( i + 1 < ${#tokens[@]} )); then
          next="${tokens[$((i + 1))]}"
        else
          next=""
        fi
        if [[ -n "$next" && ( "$next" != -* || "$next" =~ ^-[0-9] ) ]]; then
          parts+=("$(sanitize_setting_part "$key")_$(sanitize_setting_part "$next")")
          ((i+=2))
        else
          parts+=("$(sanitize_setting_part "$key")_true")
          ((i+=1))
        fi
        ;;
      *)
        # Unknown positional and short-form arguments cannot be interpreted
        # safely without knowing whether they consume another token.
        ((i+=1))
        ;;
    esac
  done

  if [[ ${#parts[@]} -eq 0 ]]; then
    printf 'default'
    return
  fi

  local result="" part
  for part in "${parts[@]}"; do
    if [[ -z "$result" ]]; then
      result="$part"
    else
      result="${result}__${part}"
    fi
  done
  printf '%s' "$result"
}

# ASTRAL-Pro3 has several meaningful short options that are unrelated to the
# STELAR-Pro short-option vocabulary above. Encode them explicitly so option
# sweeps cannot overwrite one another (for example, default rounds versus -R).
build_astral_pro3_setting_name_from_opts() {
  local raw="$1"
  local -a tokens=()
  local -a parts=()
  local i=0 token key value next

  if [[ -z "${raw// }" ]]; then
    printf 'default'
    return
  fi

  read -r -a tokens <<< "$raw"
  while (( i < ${#tokens[@]} )); do
    token="${tokens[$i]}"
    case "$token" in
      -v|--verbose)
        # Verbosity does not define an experimental setting, but it consumes a
        # value in ASTRAL-Pro3 (unlike STELAR-Pro's -v flag).
        if (( i + 1 < ${#tokens[@]} )); then ((i+=2)); else ((i+=1)); fi
        ;;
      -R|--moreround) parts+=("more-rounds_true"); ((i+=1)) ;;
      -C|--scoring) parts+=("scoring_true"); ((i+=1)) ;;
      -E|--noexit) parts+=("resolve-polytomies_true"); ((i+=1)) ;;
      -t|--thread)
        if (( i + 1 < ${#tokens[@]} )); then
          parts+=("threads_$(sanitize_setting_part "${tokens[$((i + 1))]}")")
          ((i+=2))
        else
          ((i+=1))
        fi
        ;;
      -r|--round|-s|--subsample|-u|--support|-w|--downweightrepeat|-l|--lambda|-a|--mapping|-c|--constraint|-g|--guide|-e|--exit)
        key="${token#--}"
        case "$token" in
          -r) key=round ;; -s) key=subsample ;; -u) key=support ;;
          -w) key=downweightrepeat ;; -l) key=lambda ;; -a) key=mapping ;;
          -c) key=constraint ;; -g) key=guide ;; -e) key=exit ;;
        esac
        if (( i + 1 < ${#tokens[@]} )); then
          parts+=("$(sanitize_setting_part "$key")_$(sanitize_setting_part "${tokens[$((i + 1))]}")")
          ((i+=2))
        else
          ((i+=1))
        fi
        ;;
      --*=*)
        key="${token%%=*}"
        value="${token#*=}"
        if [[ "$key" == --verbose ]]; then
          ((i+=1))
          continue
        fi
        parts+=("$(sanitize_setting_part "${key#--}")_$(sanitize_setting_part "$value")")
        ((i+=1))
        ;;
      --*)
        key="${token#--}"
        next="${tokens[$((i + 1))]:-}"
        if [[ -n "$next" && ( "$next" != -* || "$next" =~ ^-[0-9] ) ]]; then
          parts+=("$(sanitize_setting_part "$key")_$(sanitize_setting_part "$next")")
          ((i+=2))
        else
          parts+=("$(sanitize_setting_part "$key")_true")
          ((i+=1))
        fi
        ;;
      *) ((i+=1)) ;;
    esac
  done

  if [[ ${#parts[@]} -eq 0 ]]; then
    printf 'default'
    return
  fi

  local result="" part
  for part in "${parts[@]}"; do
    if [[ -z "$result" ]]; then result="$part"; else result+="__${part}"; fi
  done
  printf '%s' "$result"
}
