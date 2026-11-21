#!/usr/bin/env bash
# The pre-promotion subset-local path uses lineage scope unless overridden.
set -euo pipefail
S6_TEMP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
S6_REPO_DIR="$(cd "${S6_TEMP_DIR}/../.." && pwd -P)"
S6_ARCHIVE_DIR="${S6_REPO_DIR}/temporary-builds/s6-without-subset-local"
S6_HEAP="${STELAR_PRO_XMX:-9g}"
S6_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --xmx)
      [[ $# -ge 2 ]] || { echo '--xmx requires a heap size' >&2; exit 2; }
      S6_HEAP="$2"; shift 2 ;;
    --no-build) shift ;;
    *) S6_ARGS+=("$1"); shift ;;
  esac
done
echo 'S6 subset-local strict lineage (temporary build; S4 at <=100 species)' >&2
exec java -Xms256m "-Xmx${S6_HEAP}" \
  "-Dstelarpro.home=${S6_REPO_DIR}" \
  "-Djava.library.path=${S6_REPO_DIR}/native" \
  -cp "${S6_ARCHIVE_DIR}/stelar-pro-s6-small-s4.jar" \
  stelarx.Main "${S6_ARGS[@]}" --search-space S6 --subsetlocal
