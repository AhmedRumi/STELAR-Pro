#!/usr/bin/env bash
#
# Build (or rebuild) the machine-specific backends STELAR-Pro depends on:
#
#   ASTER-Linux/bin/astral-pro3   rooting/tagging backend  (g++ -march=native; required)
#   native/libstelar_pro_*.so     CUDA kernels             (nvcc; optional, CPU fallback)
#
# These binaries are deliberately NOT tracked by git. A binary built on one
# machine generally does not load on another: newer glibc symbol versions
# ("version `GLIBC_2.38' not found"), a different CPU (-march=native -> SIGILL)
# or a different GPU generation (compute capability below the artifact's
# minimum). run.sh runs this script before every run; it is a fast no-op when
# everything was built on this machine from the current sources.
#
# Usage: ./ensure_backends.sh [--force] [--cpu-only] [--skip-aster] [--quiet]
#   CUDA_ARCH=<arch>   nvcc target for the CUDA libraries (default: native, i.e.
#                      exactly this machine's GPU; all-major for a portable build)
#   NVCC=<path>        nvcc to use when it is not on PATH
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/backend_fingerprint.sh
source "${ROOT}/scripts/backend_fingerprint.sh"

ASTER_DIR="${ROOT}/ASTER-Linux"
ASTER_BIN="${ASTER_DIR}/bin/astral-pro3"
ASTER_STAMP="${ASTER_DIR}/bin/.build-stamp"
NATIVE_DIR="${ROOT}/native"
NATIVE_STAMP="${NATIVE_DIR}/.build-stamp"
NATIVE_LIBS=(libstelar_pro_weight.so libstelar_pro_dp.so libstelar_pro_dist.so libstelar_pro_sim.so)
PROBE_TIMEOUT_SEC=60

FORCE=false
CPU_ONLY=false
SKIP_ASTER=false
QUIET=false

usage() {
  cat <<USAGE
Usage: ./ensure_backends.sh [options]

Builds astral-pro3 and the CUDA libraries for THIS machine when they are
missing, were built on another machine, or are older than their sources.

Options:
  --force        Rebuild everything even if it looks current
  --cpu-only     Do not build the CUDA libraries
  --skip-aster   Do not manage ASTER-Linux/bin/astral-pro3 (an override is in use)
  --quiet, -q    Print only warnings and errors
  -h, --help     Show this message

Environment:
  CUDA_ARCH      nvcc architecture (default: native). Use all-major for a build
                 that must run on GPUs other than the one in this machine.
  NVCC           Path to nvcc when it is not on PATH.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --force) FORCE=true ;;
    --cpu-only|--without-cuda) CPU_ONLY=true ;;
    --skip-aster) SKIP_ASTER=true ;;
    --quiet|-q) QUIET=true ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

if [[ -t 2 && -z "${NO_COLOR:-}" ]]; then
  GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; DIM='\033[2m'; NC='\033[0m'
else
  GREEN=''; YELLOW=''; RED=''; DIM=''; NC=''
fi
say()  { [[ "$QUIET" == true ]] || echo -e "$@" >&2; }
warn() { echo -e "${YELLOW}Warning:${NC} $*" >&2; }
fail() { echo -e "${RED}Error:${NC} $*" >&2; }

HOST_GLIBC="$(stelar_glibc_version)"
HOST_CPU="$(stelar_cpu_model)"
HOST_GPUS="$(stelar_gpu_list)"

first_line() { printf '%s\n' "$1" | sed '/^[[:space:]]*$/d' | head -n1; }

with_timeout() {
  if command -v timeout >/dev/null 2>&1; then timeout "$PROBE_TIMEOUT_SEC" "$@"; else "$@"; fi
}

# Any file under $1 with one of the given suffixes newer than $2?
sources_newer_than() {
  local dir="$1" ref="$2"; shift 2
  local args=() first=true suffix
  for suffix in "$@"; do
    if [[ "$first" == true ]]; then first=false; else args+=(-o); fi
    args+=(-name "*.${suffix}")
  done
  [[ -n "$(find "$dir" -type f \( "${args[@]}" \) -newer "$ref" -print -quit 2>/dev/null)" ]]
}

# PATH for compiling against the system glibc. GCC locates `ld` through PATH,
# and an active conda/mamba environment puts its own binutils first; that linker
# cannot resolve the system libc (undefined ...@GLIBC_PRIVATE references).
toolchain_path() {
  local cleaned
  cleaned="$(printf '%s' "$PATH" | tr ':' '\n' | grep -viE 'conda|mamba' | paste -sd: -)"
  if [[ -x /usr/bin/ld && "$(PATH="$cleaned" command -v ld 2>/dev/null || true)" != /usr/bin/ld ]]; then
    cleaned="/usr/bin:/bin:${cleaned}"
  fi
  printf '%s' "$cleaned"
}

# ── ASTRAL-Pro3 (rooting/tagging backend) ────────────────────────────────────

ensure_aster() {
  if [[ "$SKIP_ASTER" == true ]]; then
    say "  ${DIM}·${NC} astral-pro3   skipped ${DIM}(--astral-pro-executable override in use)${NC}"
    return 0
  fi
  if [[ -n "${STELAR_PRO_EXECUTABLE:-}" ]]; then
    say "  ${DIM}·${NC} astral-pro3   skipped ${DIM}(STELAR_PRO_EXECUTABLE=${STELAR_PRO_EXECUTABLE})${NC}"
    return 0
  fi

  local reason="" probe_out=""
  if [[ "$FORCE" == true ]]; then
    reason="--force"
  elif [[ ! -x "$ASTER_BIN" ]]; then
    reason="not built yet"
  elif [[ ! -f "$ASTER_STAMP" ]]; then
    reason="no build record; the binary may have been built on another machine"
  elif [[ "$(stelar_stamp_get "$ASTER_STAMP" glibc)" != "$HOST_GLIBC" ]]; then
    reason="built against glibc $(stelar_stamp_get "$ASTER_STAMP" glibc), this machine has glibc ${HOST_GLIBC}"
  elif [[ "$(stelar_stamp_get "$ASTER_STAMP" cpu)" != "$HOST_CPU" ]]; then
    reason="built with -march=native for a different CPU"
  elif sources_newer_than "${ASTER_DIR}/src" "$ASTER_BIN" cpp hpp h c; then
    reason="ASTER sources are newer than the binary"
  elif ! probe_out="$(with_timeout "$ASTER_BIN" -h 2>&1)"; then
    reason="does not run here: $(first_line "$probe_out")"
  fi

  if [[ -z "$reason" ]]; then
    say "  ${GREEN}✓${NC} astral-pro3   current ${DIM}(built on this machine, glibc ${HOST_GLIBC})${NC}"
    return 0
  fi

  local build_path; build_path="$(toolchain_path)"
  if ! PATH="$build_path" command -v g++ >/dev/null 2>&1 || ! PATH="$build_path" command -v make >/dev/null 2>&1; then
    fail "astral-pro3 must be built for this machine (${reason}), but g++/make are missing."
    echo "  On Ubuntu/Debian: sudo apt install build-essential" >&2
    echo "  Or point STELAR_PRO_EXECUTABLE / --astral-pro-executable at a working ASTRAL-Pro3." >&2
    return 1
  fi

  say "  ${YELLOW}⟳${NC} astral-pro3   building for this machine: ${reason} ${DIM}(about a minute)${NC}"
  mkdir -p "${ASTER_DIR}/bin"
  local log="${ASTER_DIR}/bin/.build.log" start=$SECONDS
  if ! PATH="$build_path" make -C "$ASTER_DIR" astral-pro >"$log" 2>&1; then
    fail "building astral-pro3 failed. Last lines of ${log}:"
    tail -n 20 "$log" >&2
    if grep -q 'GLIBC_PRIVATE' "$log"; then
      echo "  A non-system linker was used (check: command -v ld). Deactivate conda (conda deactivate) or install binutils." >&2
    fi
    return 1
  fi
  if ! probe_out="$(with_timeout "$ASTER_BIN" -h 2>&1)"; then
    fail "freshly built astral-pro3 does not run: $(first_line "$probe_out")"
    return 1
  fi
  {
    echo "built_on=$(hostname) $(date -Is)"
    echo "glibc=${HOST_GLIBC}"
    echo "cpu=${HOST_CPU}"
    echo "compiler=$(g++ --version | head -n1)"
  } > "$ASTER_STAMP"
  say "  ${GREEN}✓${NC} astral-pro3   built in $((SECONDS - start)) s"
}

# ── CUDA libraries ───────────────────────────────────────────────────────────

find_nvcc() {
  if [[ -n "${NVCC:-}" && -x "${NVCC}" ]]; then echo "$NVCC"; return; fi
  if command -v nvcc >/dev/null 2>&1; then command -v nvcc; return; fi
  [[ -x /usr/local/cuda/bin/nvcc ]] && echo /usr/local/cuda/bin/nvcc
  return 0
}

native_libs_missing() {
  local lib
  for lib in "${NATIVE_LIBS[@]}"; do [[ -f "${NATIVE_DIR}/${lib}" ]] || return 0; done
  return 1
}

# Prints the first unresolved dependency of any CUDA library (empty when all load).
native_libs_unresolved() {
  local lib
  for lib in "${NATIVE_LIBS[@]}"; do
    [[ -f "${NATIVE_DIR}/${lib}" ]] || continue
    ldd "${NATIVE_DIR}/${lib}" 2>&1 | grep -m1 'not found' | sed 's/^[[:space:]]*//' && return 0
  done
  return 0
}

ensure_native() {
  if [[ "$CPU_ONLY" == true ]]; then
    say "  ${DIM}·${NC} CUDA libs     skipped ${DIM}(CPU-only run)${NC}"
    return 0
  fi
  local nvcc arch="${CUDA_ARCH:-}"
  nvcc="$(find_nvcc)"
  if [[ -z "$nvcc" ]]; then
    say "  ${DIM}·${NC} CUDA libs     not built: nvcc not found ${DIM}(CPU fallback; install the CUDA toolkit to use the GPU)${NC}"
    local unresolved; unresolved="$(native_libs_unresolved)"
    [[ -z "$unresolved" ]] || warn "native/ holds CUDA libraries from another machine that do not load here (${unresolved}); they will be ignored."
    return 0
  fi
  if [[ -z "$HOST_GPUS" && ( -z "$arch" || "$arch" == native ) ]]; then
    say "  ${DIM}·${NC} CUDA libs     not built: no NVIDIA GPU detected ${DIM}(CPU fallback; CUDA_ARCH=all-major builds anyway)${NC}"
    return 0
  fi

  local reason="" unresolved=""
  if [[ "$FORCE" == true ]]; then
    reason="--force"
  elif native_libs_missing; then
    reason="not built yet"
  elif [[ ! -f "$NATIVE_STAMP" ]]; then
    reason="no build record; the libraries may have been built on another machine"
  elif [[ "$(stelar_stamp_get "$NATIVE_STAMP" glibc)" != "$HOST_GLIBC" ]]; then
    reason="built against glibc $(stelar_stamp_get "$NATIVE_STAMP" glibc), this machine has glibc ${HOST_GLIBC}"
  elif [[ "$(stelar_stamp_get "$NATIVE_STAMP" gpus)" != "$HOST_GPUS" ]]; then
    reason="built for GPU '$(stelar_stamp_get "$NATIVE_STAMP" gpus)', this machine has '${HOST_GPUS}'"
  elif [[ -n "$arch" && "$(stelar_stamp_get "$NATIVE_STAMP" cuda_arch)" != "$arch" ]]; then
    reason="CUDA_ARCH=${arch} requested, previously built with $(stelar_stamp_get "$NATIVE_STAMP" cuda_arch)"
  elif sources_newer_than "${ROOT}/src/native" "${NATIVE_DIR}/${NATIVE_LIBS[0]}" cu cuh h hpp; then
    reason="CUDA sources are newer than the libraries"
  else
    unresolved="$(native_libs_unresolved)"
    [[ -z "$unresolved" ]] || reason="do not load here: ${unresolved}"
  fi

  if [[ -z "$reason" ]]; then
    say "  ${GREEN}✓${NC} CUDA libs     current ${DIM}(built on this machine for ${HOST_GPUS}, arch $(stelar_stamp_get "$NATIVE_STAMP" cuda_arch))${NC}"
    return 0
  fi

  say "  ${YELLOW}⟳${NC} CUDA libs     building for this machine: ${reason} ${DIM}(nvcc, may take a few minutes)${NC}"
  local log="${NATIVE_DIR}/.build.log" start=$SECONDS
  mkdir -p "$NATIVE_DIR"
  if ! CUDA_ARCH="${arch:-native}" NVCC="$nvcc" "${ROOT}/build_native.sh" >"$log" 2>&1; then
    warn "building the CUDA libraries failed; continuing with the CPU fallback. Last lines of ${log}:"
    tail -n 20 "$log" >&2
    return 0
  fi
  unresolved="$(native_libs_unresolved)"
  if [[ -n "$unresolved" ]]; then
    warn "freshly built CUDA libraries still do not load (${unresolved}); the CPU fallback will be used."
    return 0
  fi
  say "  ${GREEN}✓${NC} CUDA libs     built in $((SECONDS - start)) s ${DIM}(arch ${arch:-native})${NC}"
}

say "=== Machine-specific backends ==="
status=0
ensure_aster || status=1
ensure_native
exit "$status"
