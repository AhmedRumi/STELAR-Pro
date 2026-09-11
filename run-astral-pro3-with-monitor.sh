#!/usr/bin/env bash
# Run the bundled native ASTRAL-Pro3 on an arbitrary Newick gene-tree file and
# record the same resource/statistics sidecar used by the experiment runners.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NTFY_CHANNEL_NAME="${NTFY_CHANNEL_NAME:-anik-phylo-stx}"

INPUT_FILE=""
OUTPUT_FILE=""
ASTRAL_PRO3_BIN=""
REFERENCE_SPECIES_TREE=""
TIME_MONITOR=true
GPU_MONITOR=true
NO_NOTIFY=false
DEBUG=0
ASTRAL_PRO3_ARGS=()

print_help() {
  cat <<EOF
run-astral-pro3-with-monitor.sh - ASTRAL-Pro3 wrapper with performance monitoring

Usage: $0 --input GENE_TREES --output SPECIES_TREE [options]

Required:
  --input, -i FILE          Newick gene trees (one per line)
  --output, -o FILE         Output species tree

Optional:
  --astral-pro3-bin FILE    Executable override
                            (default: ASTER-Linux/bin/astral-pro3)
  --opts, --alg-opts "..." Extra ASTRAL-Pro3 options, such as "-t 16 --seed 42"
  --reference-species-tree FILE
                            Calculate RF rate against this tree
  --threads, -t N          Convenience alias passed to ASTRAL-Pro3
  --mapping, -a FILE       Gene-to-species mapping file
  --root TAXON             Root output at TAXON
  --no-time-monitor        Disable CPU-RAM sampling through time -v
  --no-gpu-monitor         Disable GPU-memory sampling
  --no-notify, -nn         Disable completion notification
  --debug                  Enable shell tracing
  --help, -h               Show this message

All arguments after -- are passed directly to ASTRAL-Pro3. The wrapper owns
-i/--input and -o/--output; do not repeat those in --opts.

Example:
  $0 -i data/all_gt.tre -o results/astral-pro3.tre --threads 16 --no-notify
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT_FILE="$2"; shift 2 ;;
    -o|--output) OUTPUT_FILE="$2"; shift 2 ;;
    --astral-pro3-bin|--astral-pro-bin|--executable) ASTRAL_PRO3_BIN="$2"; shift 2 ;;
    --opts|--alg-opts|--astral-pro3-opts)
      read -r -a EXTRA_ARGS <<< "$2"
      ASTRAL_PRO3_ARGS+=("${EXTRA_ARGS[@]}")
      shift 2
      ;;
    --reference-species-tree) REFERENCE_SPECIES_TREE="$2"; shift 2 ;;
    --threads|-t) ASTRAL_PRO3_ARGS+=(--thread "$2"); shift 2 ;;
    --mapping|-a) ASTRAL_PRO3_ARGS+=(--mapping "$2"); shift 2 ;;
    --root) ASTRAL_PRO3_ARGS+=(--root "$2"); shift 2 ;;
    --no-time-monitor) TIME_MONITOR=false; shift ;;
    --no-gpu-monitor) GPU_MONITOR=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --debug) DEBUG=1; shift ;;
    --help|-h) print_help; exit 0 ;;
    --) shift; ASTRAL_PRO3_ARGS+=("$@"); break ;;
    *) echo "Error: unknown option: $1" >&2; print_help >&2; exit 2 ;;
  esac
done

if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" ]]; then
  echo "Error: both --input and --output are required." >&2
  exit 2
fi

for option in "${ASTRAL_PRO3_ARGS[@]}"; do
  case "$option" in
    --search-space|--search-space=*)
      echo "Error: --search-space is a STELAR-Pro-only option; ASTRAL-Pro3 does not support it." >&2
      echo "Use --round N and --subsample N, use -R, or omit it for ASTRAL-Pro3 defaults." >&2
      exit 2
      ;;
  esac
done

INPUT_FILE="$(realpath -m "$INPUT_FILE")"
OUTPUT_FILE="$(realpath -m "$OUTPUT_FILE")"
if [[ -n "$REFERENCE_SPECIES_TREE" ]]; then
  REFERENCE_SPECIES_TREE="$(realpath -m "$REFERENCE_SPECIES_TREE")"
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: input file not found: $INPUT_FILE" >&2
  exit 3
fi
if [[ -n "$REFERENCE_SPECIES_TREE" && ! -f "$REFERENCE_SPECIES_TREE" ]]; then
  echo "Error: reference species tree not found: $REFERENCE_SPECIES_TREE" >&2
  exit 3
fi

if [[ -z "$ASTRAL_PRO3_BIN" ]]; then
  # The bundled executable is machine-specific. Keep it current without
  # needlessly building STELAR-Pro's CUDA libraries.
  "${ROOT}/ensure_backends.sh" --cpu-only --quiet
  ASTRAL_PRO3_BIN="${ROOT}/ASTER-Linux/bin/astral-pro3"
elif [[ "$ASTRAL_PRO3_BIN" != /* ]]; then
  ASTRAL_PRO3_BIN="${ROOT}/${ASTRAL_PRO3_BIN}"
fi
ASTRAL_PRO3_BIN="$(realpath -m "$ASTRAL_PRO3_BIN")"
if [[ ! -x "$ASTRAL_PRO3_BIN" ]]; then
  echo "Error: ASTRAL-Pro3 executable not found or not executable: $ASTRAL_PRO3_BIN" >&2
  exit 4
fi

[[ "$DEBUG" == 1 ]] && set -x
mkdir -p "$(dirname "$OUTPUT_FILE")"

TEMP_DIR="$(mktemp -d)"
TIME_TMP="${TEMP_DIR}/astral-pro3-time.log"
GPU_TMP="${TEMP_DIR}/astral-pro3-gpu.log"
DONE_FILE="${TEMP_DIR}/.done"
cleanup() {
  touch "$DONE_FILE" 2>/dev/null || true
  if [[ -n "${MON_PID:-}" ]]; then
    kill "$MON_PID" 2>/dev/null || true
    wait "$MON_PID" 2>/dev/null || true
  fi
  rm -rf "$TEMP_DIR" 2>/dev/null || true
}
trap cleanup EXIT

TIME_CMD=""
if [[ "$TIME_MONITOR" == true ]]; then
  if [[ -x /usr/bin/time ]]; then
    TIME_CMD=/usr/bin/time
  else
    echo "Warning: /usr/bin/time is unavailable; CPU RAM will be reported as NA." >&2
    TIME_MONITOR=false
  fi
fi

MON_PID=""
if [[ "$GPU_MONITOR" == true ]] && command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
  (
    curmax=0
    while [[ ! -f "$DONE_FILE" ]]; do
      gpu_val=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null \
        | awk 'BEGIN{m=0} {v=int($1); if(v>m) m=v} END{print m+0}')
      if [[ "$gpu_val" =~ ^[0-9]+$ ]] && (( gpu_val > curmax )); then curmax=$gpu_val; fi
      sleep 0.2
    done
    echo "$curmax" > "$GPU_TMP"
  ) &
  MON_PID=$!
else
  GPU_MONITOR=false
fi

CMD=("$ASTRAL_PRO3_BIN" "${ASTRAL_PRO3_ARGS[@]}" -i "$INPUT_FILE" -o "$OUTPUT_FILE")
printf '=== ASTRAL-Pro3 Monitor Wrapper ===\n'
printf 'Input file:       %s\n' "$INPUT_FILE"
printf 'Output file:      %s\n' "$OUTPUT_FILE"
printf 'Executable:       %s\n' "$ASTRAL_PRO3_BIN"
printf 'Options:          %s\n' "${ASTRAL_PRO3_ARGS[*]:-(defaults)}"
printf 'Time monitor:     %s\n' "$TIME_MONITOR"
printf 'GPU monitor:      %s\n\n' "$GPU_MONITOR"

START_NS=$(date +%s%N)
set +e
if [[ "$TIME_MONITOR" == true ]]; then
  "$TIME_CMD" -v -o "$TIME_TMP" "${CMD[@]}"
else
  "${CMD[@]}"
fi
ASTRAL_PRO3_EXIT_CODE=$?
set -e
touch "$DONE_FILE"
END_NS=$(date +%s%N)

if [[ -n "$MON_PID" ]]; then wait "$MON_PID" 2>/dev/null || true; fi
RUNNING_TIME=$(awk "BEGIN {printf \"%.3f\", (${END_NS}-${START_NS})/1000000000}")
MAX_CPU_MB=NA
MAX_GPU_MB=NA
if [[ -s "$TIME_TMP" ]]; then
  MAX_RSS_KB=$(awk -F: '/Maximum resident set size/ {gsub(/^[ \t]+/,"",$2); print int($2); exit}' "$TIME_TMP")
  if [[ "$MAX_RSS_KB" =~ ^[0-9]+$ ]]; then
    MAX_CPU_MB=$(awk "BEGIN {printf \"%.3f\", ${MAX_RSS_KB}/1024}")
  fi
fi
if [[ -s "$GPU_TMP" ]]; then
  MAX_GPU_VAL=$(head -n1 "$GPU_TMP")
  if [[ "$MAX_GPU_VAL" =~ ^[0-9]+$ ]]; then
    MAX_GPU_MB=$(awk "BEGIN {printf \"%.3f\", ${MAX_GPU_VAL} * 1.024}")
  fi
fi

if [[ $ASTRAL_PRO3_EXIT_CODE -eq 0 && ! -s "$OUTPUT_FILE" ]]; then
  echo "Error: ASTRAL-Pro3 exited successfully but produced no output tree." >&2
  ASTRAL_PRO3_EXIT_CODE=5
fi

RF_RATE=NA
PYTHON_BIN="${STELAR_PRO_PYTHON:-${ROOT}/.venv/bin/python}"
[[ -x "$PYTHON_BIN" ]] || PYTHON_BIN=python3
if [[ $ASTRAL_PRO3_EXIT_CODE -eq 0 && -n "$REFERENCE_SPECIES_TREE" && -f "${ROOT}/rf.py" ]]; then
  rf_output=$("$PYTHON_BIN" "${ROOT}/rf.py" "$OUTPUT_FILE" "$REFERENCE_SPECIES_TREE" 2>&1) || true
  rf_line=$(grep -i 'Robinson-Foulds distance' <<< "$rf_output" | tail -n1 || true)
  if [[ -n "$rf_line" ]]; then
    RF_RATE=$(grep -Eo '[0-9]+([.][0-9]+)?' <<< "$rf_line" | tail -n1 || echo NA)
  fi
fi

OUTPUT_BASENAME="$(basename "$OUTPUT_FILE")"
OUTPUT_STEM="$OUTPUT_BASENAME"
[[ "$OUTPUT_BASENAME" == *.* ]] && OUTPUT_STEM="${OUTPUT_BASENAME%.*}"
STATS_FILE="$(dirname "$OUTPUT_FILE")/${OUTPUT_STEM}_stats.csv"
printf '%s\n' 'algorithm,input_file,output_file,running_time_s,max_cpu_mb,max_gpu_mb,optimal_triplet_score,rf_rate,exit_code' > "$STATS_FILE"
printf 'astral-pro3,%s,%s,%s,%s,%s,NA,%s,%s\n' \
  "$(basename "$INPUT_FILE")" "$(basename "$OUTPUT_FILE")" "$RUNNING_TIME" \
  "$MAX_CPU_MB" "$MAX_GPU_MB" "$RF_RATE" "$ASTRAL_PRO3_EXIT_CODE" >> "$STATS_FILE"

echo
echo "=== ASTRAL-Pro3 Execution Summary ==="
echo "Exit code:       $ASTRAL_PRO3_EXIT_CODE"
echo "Running time:    ${RUNNING_TIME}s"
echo "Max CPU RAM:     ${MAX_CPU_MB} MB"
echo "Max GPU VRAM:    ${MAX_GPU_MB} MB"
[[ -n "$REFERENCE_SPECIES_TREE" ]] && echo "RF rate:         $RF_RATE"
echo "Stats saved to:  $STATS_FILE"

if [[ "$NO_NOTIFY" == false ]] && command -v curl >/dev/null 2>&1; then
  STATUS_EMOJI=$(if [[ $ASTRAL_PRO3_EXIT_CODE -eq 0 ]]; then echo "✅"; else echo "❌"; fi)
  STATUS_TEXT=$(if [[ $ASTRAL_PRO3_EXIT_CODE -eq 0 ]]; then echo "completed"; else echo "failed (exit $ASTRAL_PRO3_EXIT_CODE)"; fi)
  DISPLAY_OPTS="${ASTRAL_PRO3_ARGS[*]:-<default>}"
  CSV_HEADER="algorithm,input_file,output_file,running_time_s,max_cpu_mb,max_gpu_mb,optimal_triplet_score,rf_rate,exit_code"
  CSV_ROW="astral-pro3,$(basename "$INPUT_FILE"),$(basename "$OUTPUT_FILE"),${RUNNING_TIME},${MAX_CPU_MB},${MAX_GPU_MB},NA,${RF_RATE},${ASTRAL_PRO3_EXIT_CODE}"
  NOTIFY_BODY="${STATUS_EMOJI} ASTRAL-Pro3 ${STATUS_TEXT}

Input: $(basename "$INPUT_FILE")
Output: $(basename "$OUTPUT_FILE")
Options: ${DISPLAY_OPTS}
RF Rate: ${RF_RATE}
Running time: ${RUNNING_TIME}s
Max CPU RAM: ${MAX_CPU_MB} MB
Max GPU VRAM: ${MAX_GPU_MB} MB

${CSV_HEADER}
${CSV_ROW}

Stats: ${STATS_FILE}"
  curl -s -d "$NOTIFY_BODY" "https://ntfy.sh/${NTFY_CHANNEL_NAME}" >/dev/null 2>&1 || true
fi

exit "$ASTRAL_PRO3_EXIT_CODE"
