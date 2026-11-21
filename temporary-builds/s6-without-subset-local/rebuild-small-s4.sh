#!/usr/bin/env bash
# Rebuild only the archived entry point; all other pre-promotion classes stay frozen.
set -euo pipefail
S6_ARCHIVE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
S6_WORK="$(mktemp -d "${TMPDIR:-/tmp}/stelar-s6-small-s4.XXXXXX")"
trap 'rm -rf "${S6_WORK}"' EXIT
S6_ORIGINAL="${S6_ARCHIVE_DIR}/stelar-pro-s6-without-subset-local.jar"
S6_UPDATED="${S6_ARCHIVE_DIR}/stelar-pro-s6-small-s4.jar"
javac -cp "${S6_ORIGINAL}" -d "${S6_WORK}" \
  "${S6_ARCHIVE_DIR}/src/stelarx/Main.java"
cp "${S6_ORIGINAL}" "${S6_UPDATED}"
jar uf "${S6_UPDATED}" -C "${S6_WORK}" stelarx/Main.class
echo "Built ${S6_UPDATED}"
