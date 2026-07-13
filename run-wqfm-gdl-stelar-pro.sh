#!/usr/bin/env bash
# Run STELAR-Pro over WQFM-GDL estimated gene-tree folders.
#
# Expected input layout:
#   BASE_DIR/estimated_gene_trees/<setting>/<replicate>/est_g_trees0001.trees
#
# Output layout:
#   BASE_DIR/output_species_trees/<setting>/<replicate>/stelar-pro/out-stelar-pro.tree
#   BASE_DIR/output_species_trees/<setting>/<replicate>/stelar-pro/stats-stelar-pro.csv
#   BASE_DIR/output_species_trees/<setting>/<replicate>/stelar-pro/run.log

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STELAR_ROOT="$SCRIPT_DIR"
BASE_DIR="/home/aaniksahaa/research/phylogeny/datasets/phylo-simulated-datasets/wqfm-gdl-data"
INPUT_ROOT=""
OUTPUT_ROOT=""
TRUE_ROOT=""
TREE_GLOB="est_g_trees*.trees"
MODE="cpu-parallel"
EXPANSION=false
XMS="512m"
XMX="8g"
FRESH=false
DRY_RUN=false
KEEP_TEMP=false
LIMIT=""
NO_GPU_MONITOR=false
JAVA_LIB_DIR=""
LIVE_LOG=true
NO_NOTIFY=false
NTFY_TOPIC="anik-stx-pro"
NO_RF=false

declare -a TAXA_FILTERS=()
declare -a SETTING_FILTERS=()
declare -a REPLICATE_FILTERS=()

usage() {
  cat <<EOF
Usage: $0 [options]

Run STELAR-Pro on WQFM-GDL data where each replicate stores gene trees in
separate files. The script creates a temporary concatenated tree file per
replicate, runs STELAR-Pro, and writes mirrored output and stats directories.

Common options:
  --base-dir DIR        WQFM-GDL base dir (default: $BASE_DIR)
  --input-root DIR      Input root (default: BASE_DIR/estimated_gene_trees)
  --output-root DIR     Output root (default: BASE_DIR/output_species_trees)
  --true-root DIR       True tree root (default: BASE_DIR/true_trees_and_MSA)
  --taxa N              Filter settings by taxa count, e.g. 200 or 500.
                        Can be repeated or comma/semicolon separated.
  --setting NAME        Run only a setting folder, e.g. sim200_dup0.25_loss0.1.
                        Can be repeated or comma/semicolon separated.
  --replicate REP       Run only a replicate folder, e.g. 01 or 15.
                        Can be repeated or comma/semicolon separated.
  --limit N             Stop after N replicate jobs.
  --fresh               Re-run even if a successful stats file already exists.
  --dry-run             Print jobs without running STELAR-Pro.

STELAR options:
  --mode MODE           cpu, cpu-parallel, gpu, or auto (default: cpu-parallel)
  --gpu                 Same as --mode gpu.
  --cpu                 Same as --mode cpu.
  --cpu-parallel        Same as --mode cpu-parallel.
  --expansion           Enable STELAR-Pro cross-tree recombination.
  --xms SIZE            Java -Xms value (default: $XMS)
  --xmx SIZE            Java -Xmx value (default: $XMX)
  --java-lib-dir DIR    Native CUDA library dir for GPU mode
                        (default: STELAR_ROOT/cuda).
  --no-gpu-monitor      Do not sample nvidia-smi for VRAM usage.
  --no-rf               Skip RF-rate calculation against true species trees.

Input options:
  --tree-glob GLOB      File glob inside each replicate (default: $TREE_GLOB)
  --keep-temp           Keep concatenated input files under each output dir.
  --no-live-log         Save STELAR-Pro output to run.log without streaming it.
  --no-notify           Disable ntfy notification after successful replicate.
  --ntfy-topic TOPIC    ntfy.sh topic/channel (default: $NTFY_TOPIC).

Examples:
  $0 --taxa 200 --mode cpu-parallel
  $0 --setting sim200_dup0.25_loss0.1 --replicate 01 --gpu --fresh
  $0 --taxa 500 --replicate 15 --expansion --xmx 32g
EOF
}

split_filter_arg() {
  local raw="$1"
  local -n out="$2"
  local work="${raw//,/ }"
  work="${work//;/ }"
  local token
  for token in $work; do
    [[ -n "$token" ]] && out+=("$token")
  done
}

normalize_replicate() {
  local rep="$1"
  if [[ "$rep" =~ ^[0-9]+$ ]]; then
    printf "%02d" "$((10#$rep))"
  else
    printf "%s" "$rep"
  fi
}

parse_taxa_from_setting() {
  local setting="$1"
  if [[ "$setting" =~ ^sim([0-9]+)_ ]]; then
    printf "%s" "${BASH_REMATCH[1]}"
  else
    printf "NA"
  fi
}

contains_value() {
  local needle="$1"
  shift
  local value
  for value in "$@"; do
    [[ "$value" == "$needle" ]] && return 0
  done
  return 1
}

setting_selected() {
  local setting="$1"
  local taxa
  taxa="$(parse_taxa_from_setting "$setting")"

  if [[ ${#SETTING_FILTERS[@]} -gt 0 ]] && ! contains_value "$setting" "${SETTING_FILTERS[@]}"; then
    return 1
  fi
  if [[ ${#TAXA_FILTERS[@]} -gt 0 ]] && ! contains_value "$taxa" "${TAXA_FILTERS[@]}"; then
    return 1
  fi
  return 0
}

replicate_selected() {
  local rep="$1"
  if [[ ${#REPLICATE_FILTERS[@]} -eq 0 ]]; then
    return 0
  fi
  contains_value "$rep" "${REPLICATE_FILTERS[@]}"
}

join_by_space() {
  local out=""
  local item
  for item in "$@"; do
    if [[ -z "$out" ]]; then
      out="$item"
    else
      out="$out $item"
    fi
  done
  printf "%s" "$out"
}

csv_write_row() {
  local file="$1"
  shift
  local first=true
  local field
  for field in "$@"; do
    if [[ "$first" == true ]]; then
      first=false
    else
      printf "," >> "$file"
    fi
    field="${field//$'\n'/ }"
    field="${field//$'\r'/ }"
    if [[ "$field" == *","* || "$field" == *"\""* || "$field" == *" "* ]]; then
      field="${field//\"/\"\"}"
      printf '"%s"' "$field" >> "$file"
    else
      printf "%s" "$field" >> "$file"
    fi
  done
  printf "\n" >> "$file"
}

write_stats() {
  local stats_file="$1"
  shift
  printf "algorithm,dataset,setting,taxa,replicate,input_file_count,mode,expansion,rf_rate,running_time_s,stelar_time_s,max_cpu_mb,max_gpu_mb,optimal_triplet_score,normalized_triplet_score,weight_calc_ms,gpu_kernel_completed,gpu_fallback,exit_code\n" > "$stats_file"
  csv_write_row "$stats_file" "$@"
}

get_stats_value() {
  local stats_file="$1"
  local column="$2"
  awk -F, -v col="$column" '
    NR==1 {
      for (i=1; i<=NF; i++) {
        if ($i == col) c=i
      }
      next
    }
    NR==2 && c {
      print $c
      exit
    }
  ' "$stats_file"
}

send_success_notification() {
  local setting="$1"
  local replicate="$2"
  local stats_file="$3"

  if [[ "$NO_NOTIFY" == true ]]; then
    return 0
  fi
  if ! command -v curl >/dev/null 2>&1; then
    echo "WARN notification skipped: curl not found" >&2
    return 0
  fi
  if [[ ! -s "$stats_file" ]]; then
    echo "WARN notification skipped: stats file missing or empty: $stats_file" >&2
    return 0
  fi

  local csv_header csv_row runtime stelar_runtime cpu_ram gpu_vram score normalized mode expansion file_count kernel_count fallback_count rf_rate
  csv_header="$(sed -n '1p' "$stats_file")"
  csv_row="$(sed -n '2p' "$stats_file")"
  runtime="$(get_stats_value "$stats_file" "running_time_s")"
  stelar_runtime="$(get_stats_value "$stats_file" "stelar_time_s")"
  cpu_ram="$(get_stats_value "$stats_file" "max_cpu_mb")"
  gpu_vram="$(get_stats_value "$stats_file" "max_gpu_mb")"
  score="$(get_stats_value "$stats_file" "optimal_triplet_score")"
  normalized="$(get_stats_value "$stats_file" "normalized_triplet_score")"
  rf_rate="$(get_stats_value "$stats_file" "rf_rate")"
  mode="$(get_stats_value "$stats_file" "mode")"
  expansion="$(get_stats_value "$stats_file" "expansion")"
  file_count="$(get_stats_value "$stats_file" "input_file_count")"
  kernel_count="$(get_stats_value "$stats_file" "gpu_kernel_completed")"
  fallback_count="$(get_stats_value "$stats_file" "gpu_fallback")"

  local message
  message=$(
    cat <<EOF
STELAR-Pro replicate completed

Run summary
  Setting: $setting
  Replicate: $replicate
  Mode: $mode
  Expansion: $expansion
  Gene tree files: $file_count

Performance
  Total runtime: ${runtime}s
  STELAR runtime: ${stelar_runtime}s
  Peak CPU RAM: ${cpu_ram} MB
  Peak GPU VRAM: ${gpu_vram} MB

Result
  RF rate: $rf_rate
  Optimal triplet score: $score
  Normalized triplet score: $normalized
  GPU kernel completions: $kernel_count
  GPU fallbacks: $fallback_count

CSV row
$csv_header
$csv_row
EOF
  )

  curl -fsS \
    -H "Title: STELAR-Pro $setting/$replicate complete" \
    -H "Tags: stelar-pro,success" \
    -d "$message" \
    "https://ntfy.sh/$NTFY_TOPIC" >/dev/null 2>&1 || {
      echo "WARN notification failed for $setting/$replicate (topic: $NTFY_TOPIC)" >&2
      return 0
    }
}

combine_gene_tree_files() {
  local replicate_dir="$1"
  local combined_file="$2"
  shift 2
  local files=("$@")

  : > "$combined_file"
  local file
  for file in "${files[@]}"; do
    awk '
      {
        gsub(/\r/, "")
        line=$0
        sub(/^[[:space:]]+/, "", line)
        sub(/[[:space:]]+$/, "", line)
        if (line != "") {
          if (line !~ /;[[:space:]]*$/) line=line ";"
          print line
        }
      }
    ' "$file" >> "$combined_file"
  done

  if [[ ! -s "$combined_file" ]]; then
    echo "Error: concatenated input is empty for $replicate_dir" >&2
    return 1
  fi
}

gpu_library_has_pro_symbols() {
  local lib_dir="$1"
  local lib="$lib_dir/libweight_calc.so"
  [[ -f "$lib" ]] || return 1
  command -v nm >/dev/null 2>&1 || return 0
  nm -D "$lib" 2>/dev/null | grep -q "launchCompactWeightCalculationPro"
}

calculate_rf_rate() {
  local inferred_tree="$1"
  local true_tree="$2"

  if [[ "$NO_RF" == true ]]; then
    printf "NA"
    return 0
  fi
  if [[ ! -s "$inferred_tree" || ! -s "$true_tree" ]]; then
    printf "NA"
    return 0
  fi
  if [[ ! -f "$STELAR_ROOT/rf.py" ]]; then
    printf "NA"
    return 0
  fi
  if ! command -v python3 >/dev/null 2>&1; then
    printf "NA"
    return 0
  fi

  local rf_output rf_candidate
  rf_output="$(python3 "$STELAR_ROOT/rf.py" "$inferred_tree" "$true_tree" 2>&1)" || {
    echo "WARN RF calculation failed for $inferred_tree vs $true_tree" >&2
    echo "$rf_output" >&2
    printf "NA"
    return 0
  }
  rf_candidate="$(printf "%s\n" "$rf_output" | awk -F': ' '/Robinson-Foulds distance:/{print $2; exit}')"
  if [[ "$rf_candidate" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    printf "%s" "$rf_candidate"
  else
    echo "WARN could not parse RF rate from rf.py output" >&2
    echo "$rf_output" >&2
    printf "NA"
  fi
}

monitor_gpu_memory() {
  local output_file="$1"
  local done_file="$2"
  local curmax=0
  while true; do
    if command -v nvidia-smi >/dev/null 2>&1; then
      local gpu_val
      gpu_val="$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | awk 'BEGIN{m=0} {v=int($1); if(v>m) m=v} END{print m+0}')"
      if [[ "$gpu_val" =~ ^[0-9]+$ ]] && (( gpu_val > curmax )); then
        curmax="$gpu_val"
      fi
    fi
    [[ -f "$done_file" ]] && break
    sleep 0.2
  done
  printf "%s\n" "$curmax" > "$output_file"
}

run_one_job() {
  local setting="$1"
  local replicate="$2"
  local replicate_dir="$3"
  local out_dir="$4"

  local out_tree="$out_dir/out-stelar-pro.tree"
  local stats_file="$out_dir/stats-stelar-pro.csv"
  local log_file="$out_dir/run.log"
  local time_file="$out_dir/time.log"
  local true_tree="$TRUE_ROOT/$setting/$replicate/s_tree.trees"
  local taxa
  taxa="$(parse_taxa_from_setting "$setting")"

  if [[ "$FRESH" == false && -f "$out_tree" && -s "$out_tree" && -f "$stats_file" ]]; then
    local prev_exit
    prev_exit="$(awk -F, 'NR==1 {for (i=1;i<=NF;i++) if ($i=="exit_code") c=i; next} NR==2 && c {print $c}' "$stats_file" 2>/dev/null || true)"
    if [[ "$prev_exit" == "0" ]]; then
      echo "SKIP $setting/$replicate (already successful)"
      return 0
    fi
  fi

  mapfile -t tree_files < <(find "$replicate_dir" -maxdepth 1 -type f -name "$TREE_GLOB" | sort -V)
  local file_count="${#tree_files[@]}"
  if [[ "$file_count" -eq 0 ]]; then
    echo "WARN $setting/$replicate: no files matching $TREE_GLOB" >&2
    return 0
  fi

  echo "JOB $setting/$replicate: $file_count gene-tree files"
  echo "  output: $out_tree"
  echo "  stats:  $stats_file"
  echo "  log:    $log_file"
  if [[ "$DRY_RUN" == true ]]; then
    echo "  input:  $replicate_dir"
    return 0
  fi

  mkdir -p "$out_dir"
  local temp_dir
  temp_dir="$(mktemp -d)"
  local combined_input="$temp_dir/combined-gene-trees.tre"
  if [[ "$KEEP_TEMP" == true ]]; then
    combined_input="$out_dir/combined-gene-trees.tre"
  fi

  local cleanup_temp=true
  [[ "$KEEP_TEMP" == true ]] && cleanup_temp=false

  combine_gene_tree_files "$replicate_dir" "$combined_input" "${tree_files[@]}"

  local mode_arg=()
  case "$MODE" in
    cpu) mode_arg=("--cpu") ;;
    cpu-parallel) mode_arg=("--cpu-parallel") ;;
    gpu) mode_arg=("--gpu") ;;
    auto) mode_arg=() ;;
    *)
      echo "Error: invalid mode '$MODE'" >&2
      exit 1
      ;;
  esac
  local expansion_arg=()
  [[ "$EXPANSION" == true ]] && expansion_arg=("--expansion")

  local lib_dir="$JAVA_LIB_DIR"
  [[ -z "$lib_dir" ]] && lib_dir="$STELAR_ROOT/cuda"
  local jar="$STELAR_ROOT/target/stelar-pro-1.0.0-SNAPSHOT.jar"
  if [[ ! -f "$jar" ]]; then
    echo "Error: JAR not found: $jar. Run ./install.sh first." >&2
    exit 1
  fi
  if [[ "$MODE" == "gpu" ]] && ! gpu_library_has_pro_symbols "$lib_dir"; then
    echo "Error: GPU mode requested, but $lib_dir/libweight_calc.so does not export STELAR-Pro CUDA symbols." >&2
    echo "       Rebuild with ./install.sh or pass --java-lib-dir pointing to a freshly compiled CUDA library." >&2
    exit 1
  fi

  local cmd=(
    java "-Xms$XMS" "-Xmx$XMX"
    "-Dstelar.root=$STELAR_ROOT"
    "-Djava.library.path=$lib_dir"
    "-Djna.platform.library.path=$lib_dir"
    -cp "$jar"
    Main -i "$combined_input" -o "$out_tree"
  )
  cmd+=("${mode_arg[@]}")
  cmd+=("${expansion_arg[@]}")

  echo "  command: ${cmd[*]}"

  local monitor_dir="$out_dir/.monitor"
  mkdir -p "$monitor_dir"
  local done_file="$monitor_dir/done"
  local gpu_file="$monitor_dir/gpu_mib"
  rm -f "$done_file" "$gpu_file"
  local gpu_pid=""
  if [[ "$NO_GPU_MONITOR" == false && "$MODE" == "gpu" && -x "$(command -v nvidia-smi || true)" ]]; then
    monitor_gpu_memory "$gpu_file" "$done_file" &
    gpu_pid="$!"
  fi

  local start_ns end_ns exit_code
  start_ns="$(date +%s%N)"
  set +e
  if [[ "$LIVE_LOG" == true ]]; then
    if [[ -x /usr/bin/time ]]; then
      /usr/bin/time -v -o "$time_file" "${cmd[@]}" 2>&1 | tee "$log_file"
      pipe_status=("${PIPESTATUS[@]}")
      exit_code="${pipe_status[0]}"
    else
      "${cmd[@]}" 2>&1 | tee "$log_file"
      pipe_status=("${PIPESTATUS[@]}")
      exit_code="${pipe_status[0]}"
      : > "$time_file"
    fi
  else
    if [[ -x /usr/bin/time ]]; then
      /usr/bin/time -v -o "$time_file" "${cmd[@]}" > "$log_file" 2>&1
      exit_code="$?"
    else
      "${cmd[@]}" > "$log_file" 2>&1
      exit_code="$?"
      : > "$time_file"
    fi
  fi
  set -e
  end_ns="$(date +%s%N)"
  touch "$done_file"
  if [[ -n "$gpu_pid" ]]; then
    wait "$gpu_pid" 2>/dev/null || true
  fi

  local running_time_s
  running_time_s="$(awk -v s="$start_ns" -v e="$end_ns" 'BEGIN {printf "%.3f", (e-s)/1000000000}')"
  local stelar_time_s
  stelar_time_s="$(awk '/Time taken:/{v=$3} END{print v ? v : "NA"}' "$log_file")"
  local score
  score="$(awk -F': ' '/OPTIMAL_TRIPLET_SCORE:/{v=$2} END{print v ? v : "NA"}' "$log_file")"
  local normalized
  normalized="$(awk -F': ' '/NORMALIZED_TRIPLET_SCORE:/{v=$2} END{print v ? v : "NA"}' "$log_file")"
  local weight_ms
  weight_ms="$(awk '/==== MEMORY-OPTIMIZED WEIGHT CALCULATION COMPLETED ====/{flag=1; next} flag && /Processing time:/{print $3; exit}' "$log_file")"
  [[ -z "$weight_ms" ]] && weight_ms="NA"
  local max_cpu_mb="NA"
  if grep -qi "Maximum resident set size" "$time_file" 2>/dev/null; then
    max_cpu_mb="$(awk -F: '/Maximum resident set size/{gsub(/^[ \t]+/,"",$2); printf "%.3f", ($2+0)/1024; exit}' "$time_file")"
  fi
  local max_gpu_mb="NA"
  if [[ -f "$gpu_file" ]]; then
    max_gpu_mb="$(awk '{printf "%.3f", ($1+0)*1.024}' "$gpu_file")"
  fi
  local kernel_count fallback_count
  kernel_count="$(grep -c "STELAR-PRO COMPACT GPU KERNEL COMPLETED SUCCESSFULLY" "$log_file" 2>/dev/null || true)"
  fallback_count="$(grep -c "Falling back to CPU" "$log_file" 2>/dev/null || true)"
  local rf_rate="NA"
  if [[ "$exit_code" -eq 0 ]]; then
    if [[ "$NO_RF" == false ]]; then
      echo "  calculating RF against true tree: $true_tree"
    fi
    rf_rate="$(calculate_rf_rate "$out_tree" "$true_tree")"
    echo "  RF rate: $rf_rate"
  fi

  write_stats "$stats_file" \
    "stelar-pro" "wqfm-gdl" "$setting" "$taxa" "$replicate" "$file_count" \
    "$MODE" "$EXPANSION" "$rf_rate" "$running_time_s" "$stelar_time_s" \
    "$max_cpu_mb" "$max_gpu_mb" "$score" "$normalized" "$weight_ms" "$kernel_count" \
    "$fallback_count" "$exit_code"

  if [[ "$cleanup_temp" == true ]]; then
    rm -rf "$temp_dir"
  fi
  rm -rf "$monitor_dir"

  if [[ "$exit_code" -eq 0 ]]; then
    send_success_notification "$setting" "$replicate" "$stats_file"
    echo "DONE $setting/$replicate: score=$score time=${running_time_s}s stats=$stats_file"
  else
    echo "FAIL $setting/$replicate: exit=$exit_code log=$log_file stats=$stats_file" >&2
    return "$exit_code"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --input-root) INPUT_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --true-root) TRUE_ROOT="$2"; shift 2 ;;
    --taxa) split_filter_arg "$2" TAXA_FILTERS; shift 2 ;;
    --setting) split_filter_arg "$2" SETTING_FILTERS; shift 2 ;;
    --replicate)
      declare -a raw_reps=()
      split_filter_arg "$2" raw_reps
      for rep in "${raw_reps[@]}"; do
        REPLICATE_FILTERS+=("$(normalize_replicate "$rep")")
      done
      shift 2
      ;;
    --limit) LIMIT="$2"; shift 2 ;;
    --mode) MODE="${2,,}"; shift 2 ;;
    --gpu) MODE="gpu"; shift ;;
    --cpu) MODE="cpu"; shift ;;
    --cpu-parallel) MODE="cpu-parallel"; shift ;;
    --expansion|-e) EXPANSION=true; shift ;;
    --xms) XMS="$2"; shift 2 ;;
    --xmx) XMX="$2"; shift 2 ;;
    --java-lib-dir) JAVA_LIB_DIR="$2"; shift 2 ;;
    --tree-glob) TREE_GLOB="$2"; shift 2 ;;
    --fresh) FRESH=true; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    --keep-temp) KEEP_TEMP=true; shift ;;
    --no-gpu-monitor) NO_GPU_MONITOR=true; shift ;;
    --no-rf) NO_RF=true; shift ;;
    --no-live-log) LIVE_LOG=false; shift ;;
    --no-notify|-nn) NO_NOTIFY=true; shift ;;
    --ntfy-topic) NTFY_TOPIC="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "$INPUT_ROOT" ]] && INPUT_ROOT="$BASE_DIR/estimated_gene_trees"
[[ -z "$OUTPUT_ROOT" ]] && OUTPUT_ROOT="$BASE_DIR/output_species_trees"
[[ -z "$TRUE_ROOT" ]] && TRUE_ROOT="$BASE_DIR/true_trees_and_MSA"

if [[ ! -d "$INPUT_ROOT" ]]; then
  echo "Error: input root does not exist: $INPUT_ROOT" >&2
  exit 1
fi

echo "=== WQFM-GDL STELAR-Pro Runner ==="
echo "Input root:  $INPUT_ROOT"
echo "Output root: $OUTPUT_ROOT"
echo "True root:   $TRUE_ROOT"
echo "Mode:        $MODE"
echo "Expansion:   $EXPANSION"
echo "RF enabled:  $([[ "$NO_RF" == true ]] && echo false || echo true)"
echo "Tree glob:   $TREE_GLOB"
[[ ${#TAXA_FILTERS[@]} -gt 0 ]] && echo "Taxa:        $(join_by_space "${TAXA_FILTERS[@]}")"
[[ ${#SETTING_FILTERS[@]} -gt 0 ]] && echo "Settings:    $(join_by_space "${SETTING_FILTERS[@]}")"
[[ ${#REPLICATE_FILTERS[@]} -gt 0 ]] && echo "Replicates:  $(join_by_space "${REPLICATE_FILTERS[@]}")"
echo

processed=0
failed=0

mapfile -t settings < <(find "$INPUT_ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -V)
for setting in "${settings[@]}"; do
  setting_selected "$setting" || continue
  setting_dir="$INPUT_ROOT/$setting"
  mapfile -t replicates < <(find "$setting_dir" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -V)
  for replicate in "${replicates[@]}"; do
    replicate="$(normalize_replicate "$replicate")"
    replicate_selected "$replicate" || continue
    replicate_dir="$setting_dir/$replicate"
    out_dir="$OUTPUT_ROOT/$setting/$replicate/stelar-pro"
    if run_one_job "$setting" "$replicate" "$replicate_dir" "$out_dir"; then
      processed=$((processed + 1))
    else
      failed=$((failed + 1))
    fi
    if [[ -n "$LIMIT" && "$processed" -ge "$LIMIT" ]]; then
      echo "Limit reached: $LIMIT"
      echo "Processed: $processed, failed: $failed"
      exit 0
    fi
  done
done

echo
echo "=== Complete ==="
echo "Processed: $processed"
echo "Failed:    $failed"
[[ "$failed" -eq 0 ]]
