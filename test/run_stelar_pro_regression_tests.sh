#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelar-pro-tests.XXXXXX")"
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

if [[ "${STELAR_PRO_SKIP_BUILD:-0}" != 1 ]]; then
  "${ROOT}/build.sh" >/dev/null
fi
python3 "${ROOT}/test/test_stelar_pro_triplets.py"
"${ROOT}/test/test_triplet_reporting.sh"

JAVA=(java -Dstelarpro.crashDir="${WORK}/expected-failure-crash-logs" \
  -cp "${ROOT}/build" stelarx.Main)
expect_failure() {
  local label="$1"
  shift
  if "$@" >"${WORK}/${label}.log" 2>&1; then
    echo "Expected failure was accepted: ${label}" >&2
    exit 1
  fi
}

for verifier in --verify-clusters --verify-partitions --verify-dp; do
  java -cp "${ROOT}/build" stelarx.Main --cpu -q \
    -i "${ROOT}/test/input/test_5taxa.tre" "$verifier" \
    >"${WORK}/${verifier#--}.log" 2>&1
  grep -q "ALL ASSERTIONS PASSED" "${WORK}/${verifier#--}.log"
done

java -cp "${ROOT}/build" stelarx.Main --cpu -q \
  -i "${ROOT}/test/input/test_incomplete.tre" \
  -o "${WORK}/default.tre" >"${WORK}/default.log" 2>&1
test -s "${WORK}/default.tre"
grep -q "Triplet score" "${WORK}/default.log"

for preset in S2 S3; do
  expect_failure "reserved-${preset}" \
    "${JAVA[@]}" --cpu -q \
      -i "${ROOT}/test/input/test_incomplete.tre" --search-space "$preset"
  grep -qF "$preset is reserved for a future STELAR-Pro implementation" \
    "${WORK}/reserved-${preset}.log"
done

for option in --intersection-method --im --weight-intersection-method; do
  label="removed-${option#--}"
  expect_failure "$label" \
    "${JAVA[@]}" --cpu -q \
      -i "${ROOT}/test/input/test_incomplete.tre" "$option" I1
  grep -qF "$option was removed" "${WORK}/${label}.log"
done

NO_COLOR=1 "${ROOT}/stelar-pro" --no-build --version >"${WORK}/version.log" 2>&1
grep -q "STELAR-Pro  v" "${WORK}/version.log"

java -Djava.library.path="${ROOT}/native" -cp "${ROOT}/build" \
  stelarx.Main --cpu --diagnose >"${WORK}/diagnose.log" 2>&1
grep -q "weight / CUDA probe:.*loaded" "${WORK}/diagnose.log"

"${ROOT}/test/test_phylogeny_data_dir.sh"
"${ROOT}/test/test_clear_bulk_simulated.sh"
"${ROOT}/test/test_bulk_simulated_exclusions.sh"

echo "STELAR-Pro focused suite: PASS"
