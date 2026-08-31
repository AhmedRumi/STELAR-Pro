#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelar-pro-tests.XXXXXX")"
cleanup() { rm -rf "$WORK"; }
trap cleanup EXIT

echo "=== STELAR-Pro implemented-stage tests ==="
"${ROOT}/build.sh" >/dev/null

mkdir -p "${WORK}/classes"
javac -cp "${ROOT}/build" -d "${WORK}/classes" \
  "${ROOT}/test/stelarx/pro/GeneTreeRooterTaggerTest.java" \
  "${ROOT}/test/stelarx/pro/GeneTreePolytomyResolverTest.java" \
  "${ROOT}/test/stelarx/pro/DuplicateAwareCandidateTest.java" \
  "${ROOT}/test/stelarx/pro/SpeciationCandidateFilterTest.java" \
  "${ROOT}/test/stelarx/tree/GeneTreeEventTagTest.java"

CP="${ROOT}/build:${WORK}/classes"
java -cp "$CP" stelarx.pro.GeneTreeRooterTaggerTest "${WORK}/root-and-tag"
java -cp "$CP" stelarx.pro.GeneTreePolytomyResolverTest "${WORK}/polytomy-resolution"
java -cp "$CP" stelarx.pro.DuplicateAwareCandidateTest "${WORK}/duplicate-candidates"
java -cp "$CP" stelarx.tree.GeneTreeEventTagTest "${WORK}/event-tags"
java -cp "$CP" stelarx.pro.SpeciationCandidateFilterTest "${WORK}/candidate-filter"

echo "STELAR-Pro implemented-stage tests: PASS"
