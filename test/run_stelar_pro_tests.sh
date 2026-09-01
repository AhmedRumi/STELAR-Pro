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
  "${ROOT}/test/stelarx/pro/MulticopyWeightIndexTest.java" \
  "${ROOT}/test/stelarx/pro/SpeciationCandidateFilterTest.java" \
  "${ROOT}/test/stelarx/tree/GeneTreeEventTagTest.java"

CP="${ROOT}/build:${WORK}/classes"
java -cp "$CP" stelar-pro.pro.GeneTreeRooterTaggerTest "${WORK}/root-and-tag"
java -cp "$CP" stelar-pro.pro.GeneTreePolytomyResolverTest "${WORK}/polytomy-resolution"
java -cp "$CP" stelar-pro.pro.DuplicateAwareCandidateTest "${WORK}/duplicate-candidates"
java -cp "$CP" stelar-pro.pro.MulticopyWeightIndexTest "${WORK}/multicopy-weight-index"
java -cp "$CP" stelar-pro.tree.GeneTreeEventTagTest "${WORK}/event-tags"
java -cp "$CP" stelar-pro.pro.SpeciationCandidateFilterTest "${WORK}/candidate-filter"

echo "STELAR-Pro implemented-stage tests: PASS"
