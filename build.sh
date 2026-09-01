#!/bin/bash
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
SRC="$ROOT/src"
BUILD="$ROOT/build"

echo "=== Building STELAR-Pro ==="
rm -rf "$BUILD"
mkdir -p "$BUILD"

TMP_SRC_LIST="$(mktemp /tmp/stelar-pro-src.XXXXXX.txt)"
trap 'rm -f "$TMP_SRC_LIST"' EXIT
find "$SRC" -name "*.java" > "$TMP_SRC_LIST"
javac -d "$BUILD" -sourcepath "$SRC" @"$TMP_SRC_LIST"
rm -f "$TMP_SRC_LIST"
trap - EXIT

echo "Build OK -> $BUILD"
echo "Run: ./stelar-pro -i <unrooted-input.tre> -vv --no-build"
