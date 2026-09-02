#!/bin/bash
# Build all CUDA JNI shared libraries for STELAR-Pro.
#   native/libstelar_pro_weight.so  -- rooted-triplet GPU weight kernel
#   native/libstelar_pro_dp.so      -- GPU cross-tree DP transition search kernel
#   native/libstelar_pro_dist.so    -- GPU distance matrix kernel (Euler tour + RMQ)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
NVCC_BIN="${NVCC:-$(command -v nvcc || true)}"
if [[ -z "$NVCC_BIN" && -x /usr/local/cuda/bin/nvcc ]]; then
  NVCC_BIN=/usr/local/cuda/bin/nvcc
fi
if [[ -z "$NVCC_BIN" ]]; then
  echo "Error: nvcc was not found. Install CUDA or set NVCC=/path/to/nvcc." >&2
  exit 1
fi
if [[ -z "${JAVA_HOME:-}" ]]; then
  JAVAC_BIN="$(command -v javac || true)"
  if [[ -z "$JAVAC_BIN" ]]; then
    echo "Error: javac was not found. Set JAVA_HOME to a JDK 21 installation." >&2
    exit 1
  fi
  JAVAC_BIN="$(readlink -f "$JAVAC_BIN")"
  JAVA_HOME="$(cd "$(dirname "$JAVAC_BIN")/.." && pwd)"
fi

# Portable default: native code for every major GPU generation supported by the
# installed toolkit, plus forward-compatible PTX for its newest generation.
# Override for a faster developer build, e.g. CUDA_ARCH=native or sm_86.
CUDA_ARCH="${CUDA_ARCH:-all-major}"
[[ "$CUDA_ARCH" == "portable" ]] && CUDA_ARCH="all-major"

# Embed the oldest compute capability actually covered by this artifact.  The
# runtime probe uses it to reject an older GPU before any kernel is launched,
# allowing Java to select the exact CPU path instead of failing mid-analysis.
MIN_CUDA_CC=0
if [[ "$CUDA_ARCH" == "all-major" ]]; then
  MIN_CUDA_CC="$("$NVCC_BIN" --list-gpu-arch | sed -n 's/^compute_\([0-9][0-9]*\)$/\1/p' | sort -n | head -1)"
elif [[ "$CUDA_ARCH" =~ ^(sm|compute)_([0-9]+)$ ]]; then
  MIN_CUDA_CC="${BASH_REMATCH[2]}"
fi
[[ -n "$MIN_CUDA_CC" ]] || MIN_CUDA_CC=0

NATIVE_OUT_DIR="${NATIVE_OUT_DIR:-${ROOT}/native}"

mkdir -p "$NATIVE_OUT_DIR"

NVCC_FLAGS=(
  -arch="${CUDA_ARCH}"
  -O3
  -Xcompiler '-fPIC'
  --shared
  -I"${JAVA_HOME}/include"
  -I"${JAVA_HOME}/include/linux"
  -DSTELAR_PRO_MIN_CUDA_CC="${MIN_CUDA_CC}"
)

# nvcc links cudart statically by default. Also remove target-machine
# libstdc++/libgcc requirements when the host toolchain supports it; glibc and
# the NVIDIA driver remain the only normal Linux runtime dependencies.
if [[ "$(uname -s)" == "Linux" ]]; then
  NVCC_FLAGS+=( -Xcompiler=-static-libstdc++ -Xcompiler=-static-libgcc )

  # CUDA 13.1 declares rsqrt/rsqrtf before glibc 2.43's C23 declarations,
  # whose added noexcept specifier otherwise makes nvcc reject every .cu file.
  # The STELAR-Pro sources do not require GNU-only libc declarations.
  MULTIARCH="$(gcc -print-multiarch 2>/dev/null || true)"
  MATHCALLS="/usr/include/${MULTIARCH}/bits/mathcalls.h"
  if "$NVCC_BIN" --version | grep -q 'release 13\.1' \
      && [[ -f "$MATHCALLS" ]] \
      && grep -q '__MATHCALL_VEC (rsqrt' "$MATHCALLS"; then
    NVCC_FLAGS+=( -U_GNU_SOURCE )
  fi
fi

echo "=== Building STELAR-Pro native GPU libraries ==="
echo "  nvcc        : $NVCC_BIN"
echo "  JDK         : $JAVA_HOME"
echo "  CUDA arch   : $CUDA_ARCH"
echo "  Minimum CC  : $MIN_CUDA_CC"
echo "  Output      : $NATIVE_OUT_DIR"

# ── Weight kernel ─────────────────────────────────────────────────────────────
SRC_W="$ROOT/src/native/stelarx_weight.cu"
OUT_W="$NATIVE_OUT_DIR/libstelar_pro_weight.so"
echo "  Building    : $SRC_W  ->  $OUT_W"
"$NVCC_BIN" "${NVCC_FLAGS[@]}" -o "$OUT_W" "$SRC_W"
echo "  OK"

# ── DP cross-tree search kernel ───────────────────────────────────────────────
SRC_DP="$ROOT/src/native/stelarx_dp.cu"
OUT_DP="$NATIVE_OUT_DIR/libstelar_pro_dp.so"
echo "  Building    : $SRC_DP  ->  $OUT_DP"
"$NVCC_BIN" "${NVCC_FLAGS[@]}" -o "$OUT_DP" "$SRC_DP"
echo "  OK"

# ── Distance matrix kernel ────────────────────────────────────────────────────
SRC_DM="$ROOT/src/native/stelarx_dist.cu"
OUT_DM="$NATIVE_OUT_DIR/libstelar_pro_dist.so"
echo "  Building    : $SRC_DM  ->  $OUT_DM"
"$NVCC_BIN" "${NVCC_FLAGS[@]}" -o "$OUT_DM" "$SRC_DM"
echo "  OK"

# ── Similarity matrix kernel ──────────────────────────────────────────────────
SRC_SIM="$ROOT/src/native/stelarx_similarity.cu"
OUT_SIM="$NATIVE_OUT_DIR/libstelar_pro_sim.so"
echo "  Building    : $SRC_SIM  ->  $OUT_SIM"
"$NVCC_BIN" "${NVCC_FLAGS[@]}" -o "$OUT_SIM" "$SRC_SIM"
echo "  OK"

echo "=== Native build complete ==="
echo "Run with:"
echo "  ./run.sh -i <input.tre> -o <output.tre> --gpu-strict -vv --no-build"
