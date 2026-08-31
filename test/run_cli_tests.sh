#!/usr/bin/env bash
# Current STELAR-X command/identity smoke tests. The superseded quartet-era
# command suite is retained under raw-prev for migration archaeology only.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-cli-tests.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

"${ROOT}/build.sh" >/dev/null

[[ -f "${ROOT}/build/stelarx/Main.class" ]]
[[ "$(find "${ROOT}/build" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')" == "stelarx" ]]

VERSION_TEXT="$(NO_COLOR=1 "${ROOT}/stelarx" --version --no-build)"
[[ "$VERSION_TEXT" == *"STELAR-X  v1.0.0"* ]]
[[ "$VERSION_TEXT" == *"Welcome to STELAR-X version 1.0.0!"* ]]

HELP_TEXT="$(NO_COLOR=1 "${ROOT}/stelarx" --help 2>&1)"
[[ "$HELP_TEXT" == *"STELAR-X wrapper"* ]]
[[ "$HELP_TEXT" == *"--search-space"* ]]
[[ "$HELP_TEXT" == *"--intersection-method"* ]]
[[ "$HELP_TEXT" == *"--gpu-strict"* ]]
[[ "$HELP_TEXT" == *"--tag-only"* ]]

DIAG_TEXT="$(NO_COLOR=1 java -Djava.library.path="${ROOT}/native" \
  -cp "${ROOT}/build" stelarx.Main --cpu --diagnose 2>&1)"
[[ "$DIAG_TEXT" == *"STELAR-X DIAGNOSTICS"* ]]
[[ "$DIAG_TEXT" == *"Selected compute:            CPU"* ]]

OUTPUT_TREE="${WORK}/species-tree.tre"
STELARX_CRASH_DIR="${WORK}/launcher-crash_logs" \
  NO_COLOR=1 "${ROOT}/stelarx" --no-build --cpu -q \
  -i "${ROOT}/test/input/stelar_candidate_5taxa.tre" \
  -o "$OUTPUT_TREE" --search-space S1 --intersection-method I2 >/dev/null
[[ -s "$OUTPUT_TREE" ]]
[[ -d "${WORK}/launcher-crash_logs" ]]

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
TAG_TEXT="$(NO_COLOR=1 "${ROOT}/stelarx" --no-build -T -q \
  -i "$TAG_INPUT" -o "$TAG_OUTPUT" \
  --astral-pro-executable "$FAKE_ASTRAL" 2>&1)"
[[ "$TAG_TEXT" == *"STELAR-Pro tag-only: rooting and tagging gene trees..."* ]]
[[ "$TAG_TEXT" == *"STELAR-Pro tag-only: wrote 1 rooted/tagged gene tree(s)"* ]]
[[ "$TAG_TEXT" != *"ASTRAL-PRO BACKEND NOISE"* ]]
grep -q ')D' "$TAG_OUTPUT"

STELARX_SKIP_BUILD=1 "${ROOT}/test/run_stelarx_tests.sh" >/dev/null

echo "STELAR-X CLI and identity tests: PASS"
