#!/usr/bin/env bash
# Current STELAR-Pro command/identity smoke tests. The superseded quartet-era
# command suite is retained under raw-prev for migration archaeology only.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelar-pro-cli-tests.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

"${ROOT}/build.sh" >/dev/null
"${ROOT}/ensure_backends.sh" --cpu-only --quiet

[[ -f "${ROOT}/build/stelarx/Main.class" ]]
[[ "$(find "${ROOT}/build" -mindepth 1 -maxdepth 1 -type d -printf '%f\n')" == "stelarx" ]]

VERSION_TEXT="$(NO_COLOR=1 "${ROOT}/stelar-pro" --version --no-build)"
[[ "$VERSION_TEXT" == *"STELAR-Pro  v1.0.0"* ]]
[[ "$VERSION_TEXT" == *"Welcome to STELAR-Pro version 1.0.0!"* ]]

HELP_TEXT="$(NO_COLOR=1 "${ROOT}/stelar-pro" --help 2>&1)"
[[ "$HELP_TEXT" == *"STELAR-Pro wrapper"* ]]
[[ "$HELP_TEXT" == *"--search-space"* ]]
[[ "$HELP_TEXT" == *"S2 (lineage-scoped DISCO/UPGMA candidates)"* ]]
[[ "$HELP_TEXT" == *"S3 (global DISCO + S1 + greedy consensus)"* ]]
[[ "$HELP_TEXT" == *"S4 (S2 candidates with rooted-triplet similarity)"* ]]
[[ "$HELP_TEXT" == *"S5 (A-Pro greedy trees + S1; rooted-triplet DP)"* ]]
[[ "$HELP_TEXT" == *"--exhaustive"* ]]
[[ "$HELP_TEXT" != *"--full-disco-candidates"* ]]
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

S2_TREE="${WORK}/s2-species-tree.tre"
S2_LOG="${WORK}/s2.log"
NO_COLOR=1 "${ROOT}/stelar-pro" --no-build --cpu \
  -i "${ROOT}/test/input/stelar_candidate_5taxa.tre" \
  -o "$S2_TREE" --search-space S2 >"$S2_LOG" 2>&1
[[ -s "$S2_TREE" ]]
grep -Fq 'S2 search: S1 candidates excluded' "$S2_LOG"
S3_TREE="${WORK}/s3-species-tree.tre"
NO_COLOR=1 "${ROOT}/stelar-pro" --no-build --cpu -q \
  -i "${ROOT}/test/input/stelar_candidate_5taxa.tre" \
  -o "$S3_TREE" --search-space S3 >/dev/null 2>&1
[[ -s "$S3_TREE" ]]
S4_TREE="${WORK}/s4-species-tree.tre"
S4_LOG="${WORK}/s4.log"
NO_COLOR=1 "${ROOT}/stelar-pro" --no-build --cpu \
  -i "${ROOT}/test/input/stelar_candidate_5taxa.tre" \
  -o "$S4_TREE" --search-space S4 >"$S4_LOG" 2>&1
[[ -s "$S4_TREE" ]]
[[ "$(grep -Fc 'S4 similarity: one rooted-triplet matrix' "$S4_LOG")" -eq 1 ]]
grep -Fq 'S4 streamed candidates:' "$S4_LOG"
grep -Fq 'S4 search: S1 candidates excluded' "$S4_LOG"
! grep -Fq 'Phase 3.5 Greedy consensus' "$S4_LOG"
! grep -Fq 'Phase 5b Cross-tree transitions' "$S4_LOG"
! grep -Fq 'Building similarity matrix (CPU)' "$S4_LOG"
for mode in S2 S4; do
  exhaustive_log="${WORK}/${mode}-exhaustive.log"
  exhaustive_tree="${WORK}/${mode}-exhaustive.tre"
  NO_COLOR=1 "${ROOT}/stelar-pro" --no-build --cpu \
    -i "${ROOT}/test/input/stelar_candidate_5taxa.tre" \
    -o "$exhaustive_tree" --search-space "$mode" --exhaustive \
    >"$exhaustive_log" 2>&1
  [[ -s "$exhaustive_tree" ]]
  grep -Fq 'all species + S1 candidates (--exhaustive)' "$exhaustive_log"
  ! grep -Fq 'S1 candidates excluded' "$exhaustive_log"
  ! grep -Fq 'Phase 5b Cross-tree transitions' "$exhaustive_log"
done
expect_launcher_failure exhaustive-s1 '--exhaustive applies only to S2 and S4' \
  --exhaustive
expect_launcher_failure removed-full-disco 'Unknown option' \
  --full-disco-candidates
S4_CROSS_LOG="${WORK}/s4-cross.log"
NO_COLOR=1 "${ROOT}/stelar-pro" --no-build --cpu \
  -i "${ROOT}/test/input/stelar_candidate_5taxa.tre" \
  -o "${WORK}/s4-cross.tre" --search-space S4 --cross-tree-splits \
  >"$S4_CROSS_LOG" 2>&1
[[ -s "${WORK}/s4-cross.tre" ]]
grep -Fq 'Phase 5b Cross-tree transitions' "$S4_CROSS_LOG"
javac -cp "${ROOT}/build" -d "$WORK" \
  "${ROOT}/test/TripletSimilarityTest.java"
java -cp "${ROOT}/build:${WORK}" TripletSimilarityTest >/dev/null 2>&1
S5_TREE="${WORK}/s5-species-tree.tre"
NO_COLOR=1 "${ROOT}/stelar-pro" --no-build --cpu -q \
  -i "${ROOT}/test/input/stelar_candidate_5taxa.tre" \
  -o "$S5_TREE" --search-space S5 --s5-rounds 2 --s5-seed 233 \
  >/dev/null 2>&1
[[ -s "$S5_TREE" ]]
javac -cp "${ROOT}/build" -d "$WORK" \
  "${ROOT}/test/stelarx/pro/S5RootInvariantCandidateTest.java"
java -cp "${ROOT}/build:${WORK}" stelarx.pro.S5RootInvariantCandidateTest \
  >/dev/null 2>&1
javac -cp "${ROOT}/build" -d "$WORK" \
  "${ROOT}/test/AProGreedySearchTest.java"
java -Djava.library.path="${ROOT}/native" -cp "${ROOT}/build:${WORK}" \
  AProGreedySearchTest >/dev/null 2>&1
java -Djava.library.path="${ROOT}/native" -cp "${ROOT}/build:${WORK}" \
  AProGreedySearchTest --large-smoke >/dev/null 2>&1
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
  '[ "$1" != "-h" ] || exit 0' \
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
