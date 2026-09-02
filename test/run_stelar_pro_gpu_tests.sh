#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelar-pro-gpu-tests.XXXXXX")"
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

if [[ "${STELAR_PRO_SKIP_BUILD:-0}" != 1 ]]; then
  "${ROOT}/build.sh" >/dev/null
fi
if [[ ! -f "${ROOT}/native/libstelar_pro_weight.so" ]]; then
  "${ROOT}/build_native.sh" >/dev/null
fi

JAVA=(java -Djava.library.path="${ROOT}/native" -cp "${ROOT}/build" stelarx.Main)

"${JAVA[@]}" --gpu-strict --diagnose >"${WORK}/diagnose.log" 2>&1
grep -q "CUDA.*" "${WORK}/diagnose.log"
grep -q "Usable:.*yes" "${WORK}/diagnose.log"

score_gpu() {
  local genes="$1" species="$2" expected="$3" log="$4"
  shift 4
  "${JAVA[@]}" --gpu-strict -q -i "$genes" \
    --score-species-tree "$species" \
    "$@" \
    >"$log" 2>&1
  grep -q "\[STELAR-Pro GPU\] weight" "$log"
  [[ "$(sed -n 's/^TRIPLET_SCORE: //p' "$log" | tail -1)" == "$expected" ]]
}

score_gpu "${ROOT}/test/input/test_5taxa.tre" \
  "${ROOT}/test/input/stelar_candidate_5taxa.tre" 21 "${WORK}/small.log"
score_gpu "${ROOT}/test/input/stelar_polytomy_5taxa.tre" \
  "${ROOT}/test/input/stelar_candidate_5taxa.tre" 11 "${WORK}/polytomy.log"
score_gpu "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre" \
  "${ROOT}/test/input/stelar_candidate_6taxa.tre" 22 \
  "${WORK}/incomplete-polytomy.log"
score_gpu "${ROOT}/all_gt_bs_rooted_37.tre" "${ROOT}/true_37.tre" \
  1390544 "${WORK}/large.log"

# The built-in scorer must retain exact results under all batching controls.
score_gpu "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre" \
  "${ROOT}/test/input/stelar_candidate_6taxa.tre" 22 \
  "${WORK}/batch-size.log" --gpu-batch-size 1
score_gpu "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre" \
  "${ROOT}/test/input/stelar_candidate_6taxa.tre" 22 \
  "${WORK}/batch-count.log" --gpu-batches 2
score_gpu "${ROOT}/test/input/stelar_polytomy_incomplete_6taxa.tre" \
  "${ROOT}/test/input/stelar_candidate_6taxa.tre" 22 \
  "${WORK}/batch-off.log" --no-gpu-batch

# Default-S1 end-to-end inference must be byte-identical on CPU and strict GPU.
"${JAVA[@]}" --cpu -q -i "${ROOT}/test/input/test_incomplete.tre" \
  -o "${WORK}/default-cpu.tre" >"${WORK}/default-cpu.log" 2>&1
"${JAVA[@]}" --gpu-strict -q -i "${ROOT}/test/input/test_incomplete.tre" \
  -o "${WORK}/default-gpu.tre" >"${WORK}/default-gpu.log" 2>&1
cmp "${WORK}/default-cpu.tre" "${WORK}/default-gpu.tre"
grep -q "Phase 6  Weight calculation.*\[GPU\]" "${WORK}/default-gpu.log"
cpu_score="$(sed -n 's/.*Triplet score[[:space:]]*\([0-9][0-9]*\).*/\1/p' \
  "${WORK}/default-cpu.log" | tail -1)"
gpu_score="$(sed -n 's/.*Triplet score[[:space:]]*\([0-9][0-9]*\).*/\1/p' \
  "${WORK}/default-gpu.log" | tail -1)"
[[ "$cpu_score" == 122 ]]
[[ "$gpu_score" == "$cpu_score" ]]

# Seeded independent oracles exercise many additional binary, incomplete, and
# internally polytomous layouts, including forced DOUBLE and INT128 kernels.
python3 "${ROOT}/test/test_stelar_pro_differential.py" \
  --gpu --cases "${STELAR_PRO_GPU_RANDOM_CASES:-6}" --no-build
python3 "${ROOT}/test/test_stelar_pro_inference.py" --gpu --no-build

echo "STELAR-Pro strict CUDA suite: PASS"
