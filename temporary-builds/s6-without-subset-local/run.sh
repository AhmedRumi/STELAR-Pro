#!/usr/bin/env bash
# Pre-promotion S6 with the small-species S4 switch; the original jar is retained.
set -euo pipefail
S6_ARCHIVE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
S6_REPO_DIR="$(cd "${S6_ARCHIVE_DIR}/../.." && pwd -P)"
S6_HEAP="${STELAR_PRO_XMX:-9g}"
S6_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --xmx)
      [[ $# -ge 2 ]] || { echo '--xmx requires a heap size' >&2; exit 2; }
      S6_HEAP="$2"; shift 2 ;;
    --no-build) shift ;;
    --subsetlocal)
      echo 'This archive is S6 without subset-local. Use production S6 for subset-local search.' >&2
      exit 2 ;;
    *) S6_ARGS+=("$1"); shift ;;
  esac
done
echo 'S6 without subset-local (temporary build; S4 at <=100 species)' >&2
exec java -Xms256m "-Xmx${S6_HEAP}" \
  "-Dstelarpro.home=${S6_REPO_DIR}" \
  "-Djava.library.path=${S6_REPO_DIR}/native" \
  -cp "${S6_ARCHIVE_DIR}/stelar-pro-s6-small-s4.jar" \
  stelarx.Main "${S6_ARGS[@]}" --search-space S6
