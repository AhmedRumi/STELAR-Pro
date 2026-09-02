#!/usr/bin/env bash
# Current STELAR-Pro command/identity smoke tests. The superseded quartet-era
# command suite is retained under raw-prev for migration archaeology only.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelar-pro-cli-tests.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

"${ROOT}/build.sh" >/dev/null

[[ -f "${ROOT}/build/stelarx/Main.class" ]]
[[ "$(find "${ROOT}/build" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')" == "stelarx" ]]

VERSION_TEXT="$(NO_COLOR=1 "${ROOT}/stelar-pro" --version --no-build)"
[[ "$VERSION_TEXT" == *"STELAR-Pro  v1.0.0"* ]]
[[ "$VERSION_TEXT" == *"Welcome to STELAR-Pro version 1.0.0!"* ]]

HELP_TEXT="$(NO_COLOR=1 "${ROOT}/stelar-pro" --help 2>&1)"
[[ "$HELP_TEXT" == *"STELAR-Pro wrapper"* ]]
[[ "$HELP_TEXT" == *"--search-space"* ]]
[[ "$HELP_TEXT" == *"S2/S3 are reserved"* ]]
[[ "$HELP_TEXT" != *"--intersection-method"* ]]
[[ "$HELP_TEXT" != *"--weight-intersection-method"* ]]
[[ "$HELP_TEXT" == *"--gpu-strict"* ]]
[[ "$HELP_TEXT" == *"--tag-only"* ]]

DIAG_TEXT="$(NO_COLOR=1 java -Djava.library.path="${ROOT}/native" \
  -cp "${ROOT}/build" stelarx.Main --cpu --diagnose 2>&1)"
[[ "$DIAG_TEXT" == *"STELAR-Pro DIAGNOSTICS"* ]]
[[ "$DIAG_TEXT" == *"Selected compute:            CPU"* ]]

OUTPUT_TREE="${WORK}/species-tree.tre"
STELAR_PRO_CRASH_DIR="${WORK}/launcher-crash_logs" \
  NO_COLOR=1 "${ROOT}/stelar-pro" --no-build --cpu -q \
  -i "${ROOT}/test/input/stelar_candidate_5taxa.tre" \
  -o "$OUTPUT_TREE" >/dev/null
[[ -s "$OUTPUT_TREE" ]]
[[ -d "${WORK}/launcher-crash_logs" ]]

expect_launcher_failure() {
  local label="$1" expected="$2"
  shift 2
  if env NO_COLOR=1 STELAR_PRO_CRASH_DIR="${WORK}/launcher-crash_logs" \
      "${ROOT}/stelar-pro" --no-build --cpu -q \
      -i "${ROOT}/test/input/stelar_candidate_5taxa.tre" \
      "$@" >"${WORK}/reject-${label}.log" 2>&1; then
    echo "Expected launcher failure was accepted: ${label}" >&2
    exit 1
  fi
  grep -Fq -- "$expected" "${WORK}/reject-${label}.log"
}

for preset in S2 S3; do
  expect_launcher_failure "reserved-${preset}" \
    "${preset} is reserved for a future STELAR-Pro implementation" \
    --search-space "$preset"
done
for option in --intersection-method --im --weight-intersection-method; do
  expect_launcher_failure "removed-${option#--}" \
    "${option} was removed" "$option" I1
done

TAG_INPUT="${WORK}/unrooted-multicopy.tre"
TAG_OUTPUT="${WORK}/rooted-tagged.tre"
FAKE_ASTRAL="${WORK}/fake-astral-pro3"
printf '(A,A,(B,C));\n' >"$TAG_INPUT"
printf '%s\n' \
  '#!/bin/sh' \
  'set -eu' \
  'echo "ASTRAL-PRO BACKEND NOISE"' \
  '[ "$1" = "-T" ]' \
  '[ "$2" = "-o" ]' \
  'printf "(B,(C,(A,A)D));\\n" > "$3"' \
  >"$FAKE_ASTRAL"
chmod +x "$FAKE_ASTRAL"
TAG_TEXT="$(NO_COLOR=1 "${ROOT}/stelar-pro" --no-build -T -q \
  -i "$TAG_INPUT" -o "$TAG_OUTPUT" \
  --astral-pro-executable "$FAKE_ASTRAL" 2>&1)"
[[ "$TAG_TEXT" == *"STELAR-Pro tag-only: rooting and tagging gene trees..."* ]]
[[ "$TAG_TEXT" == *"STELAR-Pro tag-only: wrote 1 rooted/tagged gene tree(s)"* ]]
[[ "$TAG_TEXT" != *"ASTRAL-PRO BACKEND NOISE"* ]]
grep -q ')D' "$TAG_OUTPUT"

STELAR_PRO_SKIP_BUILD=1 "${ROOT}/test/run_stelar_pro_tests.sh" >/dev/null

echo "STELAR-Pro CLI and identity tests: PASS"
