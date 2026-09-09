#!/usr/bin/env bash
# Shared helpers describing what a machine-specific backend binary depends on.
# Sourced by build_native.sh (to record a build stamp) and ensure_backends.sh
# (to decide whether an existing binary was built for this machine).

# glibc version, e.g. 2.39. Symbol versions are the usual reason a binary
# built on a newer distribution refuses to load on an older one.
stelar_glibc_version() {
  local v
  v="$(getconf GNU_LIBC_VERSION 2>/dev/null | awk '{print $2}')" || true
  if [[ -z "$v" ]]; then
    v="$(ldd --version 2>/dev/null | head -n1 | grep -oE '[0-9]+\.[0-9]+' | tail -n1)" || true
  fi
  echo "${v:-unknown}"
}

# CPU model string; ASTER is compiled with -march=native, so a binary built on
# a different CPU can die with SIGILL.
stelar_cpu_model() {
  local m
  m="$(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ *//')" || true
  echo "${m:-unknown}"
}

# Sorted "name cc X.Y" list of visible NVIDIA GPUs (empty when none/driver absent).
stelar_gpu_list() {
  command -v nvidia-smi >/dev/null 2>&1 || return 0
  nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader 2>/dev/null \
    | sed 's/, */ cc /' | sort -u | paste -sd';' - || true
}

# Read one key from a key=value stamp file (empty when absent).
stelar_stamp_get() {
  local file="$1" key="$2"
  [[ -f "$file" ]] || return 0
  sed -n "s/^${key}=//p" "$file" | head -n1
}
