#!/usr/bin/env bash
# Collect WQFM-GDL STELAR-Pro per-replicate stats into merged and summary CSVs.

set -euo pipefail

BASE_DIR="/home/aaniksahaa/research/phylogeny/datasets/phylo-simulated-datasets/wqfm-gdl-data"
OUTPUT_ROOT=""
MERGED_CSV=""
SUMMARY_CSV=""

declare -a TAXA_FILTERS=()
declare -a SETTING_FILTERS=()
declare -a REPLICATE_FILTERS=()

usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  --base-dir DIR        WQFM-GDL base dir (default: $BASE_DIR)
  --output-root DIR     Output root (default: BASE_DIR/output_species_trees)
  --output FILE         Merged CSV path
                        (default: OUTPUT_ROOT/stelar-pro-wqfm-gdl-all.csv)
  --summary-output FILE Summary CSV path
                        (default: OUTPUT_ROOT/stelar-pro-wqfm-gdl-summary.csv)
  --taxa N              Keep only taxa count N. Repeat or comma/semicolon separate.
  --setting NAME        Keep only one setting. Repeat or comma/semicolon separate.
  --replicate REP       Keep only one replicate, e.g. 01 or 15.
  --help, -h            Show this help.

The merged CSV contains one row per replicate. The summary CSV aggregates by
algorithm, dataset, setting, taxa, mode, and expansion.
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

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-dir|-b) BASE_DIR="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --output|-o) MERGED_CSV="$2"; shift 2 ;;
    --summary-output) SUMMARY_CSV="$2"; shift 2 ;;
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
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -z "$OUTPUT_ROOT" ]] && OUTPUT_ROOT="$BASE_DIR/output_species_trees"
[[ -z "$MERGED_CSV" ]] && MERGED_CSV="$OUTPUT_ROOT/stelar-pro-wqfm-gdl-all.csv"
[[ -z "$SUMMARY_CSV" ]] && SUMMARY_CSV="$OUTPUT_ROOT/stelar-pro-wqfm-gdl-summary.csv"

if [[ ! -d "$OUTPUT_ROOT" ]]; then
  echo "Error: output root does not exist: $OUTPUT_ROOT" >&2
  exit 1
fi

mkdir -p "$(dirname "$MERGED_CSV")" "$(dirname "$SUMMARY_CSV")"

filter_csv_rows() {
  local file="$1"
  awk -F, \
    -v taxa_list=" ${TAXA_FILTERS[*]} " \
    -v setting_list=" ${SETTING_FILTERS[*]} " \
    -v replicate_list=" ${REPLICATE_FILTERS[*]} " '
    NR==1 {
      for (i=1; i<=NF; i++) {
        idx[$i]=i
      }
      next
    }
    NR==2 {
      taxa=$idx["taxa"]
      setting=$idx["setting"]
      replicate=$idx["replicate"]
      if (taxa_list != "  " && index(taxa_list, " " taxa " ") == 0) next
      if (setting_list != "  " && index(setting_list, " " setting " ") == 0) next
      if (replicate_list != "  " && index(replicate_list, " " replicate " ") == 0) next
      print
    }
  ' "$file"
}

mapfile -t stat_files < <(find "$OUTPUT_ROOT" -path "*/stelar-pro/stats-stelar-pro.csv" -type f | sort -V)
if [[ "${#stat_files[@]}" -eq 0 ]]; then
  echo "Error: no stats-stelar-pro.csv files found under $OUTPUT_ROOT" >&2
  exit 1
fi

header=""
rows=0
: > "$MERGED_CSV"
for stat_file in "${stat_files[@]}"; do
  if [[ ! -s "$stat_file" ]]; then
    echo "WARN skipping empty stats file: $stat_file" >&2
    continue
  fi
  current_header="$(head -n 1 "$stat_file")"
  if [[ -z "$header" ]]; then
    header="$current_header"
    printf "%s\n" "$header" > "$MERGED_CSV"
  elif [[ "$current_header" != "$header" ]]; then
    echo "WARN skipping stats file with different header: $stat_file" >&2
    continue
  fi

  while IFS= read -r row; do
    [[ -z "$row" ]] && continue
    printf "%s\n" "$row" >> "$MERGED_CSV"
    rows=$((rows + 1))
  done < <(filter_csv_rows "$stat_file")
done

if [[ "$rows" -eq 0 ]]; then
  echo "Error: no rows matched the requested filters" >&2
  exit 1
fi

awk -F, '
  NR==1 {
    for (i=1; i<=NF; i++) idx[$i]=i
    print "algorithm,dataset,setting,taxa,mode,expansion,total_replicates,successful_replicates,failed_replicates,avg_rf_rate,avg_running_time_s,avg_stelar_time_s,avg_max_cpu_mb,avg_max_gpu_mb,avg_weight_calc_ms,avg_optimal_triplet_score"
    next
  }
  NR>1 {
    key=$idx["algorithm"] "," $idx["dataset"] "," $idx["setting"] "," $idx["taxa"] "," $idx["mode"] "," $idx["expansion"]
    total[key]++
    if ($idx["exit_code"] == 0) success[key]++; else failed[key]++
    run[key]+=$idx["running_time_s"]+0
    stelar[key]+=$idx["stelar_time_s"]+0
    cpu[key]+=$idx["max_cpu_mb"]+0
    if ($idx["max_gpu_mb"] != "NA") {
      gpu[key]+=$idx["max_gpu_mb"]+0
      gpu_count[key]++
    }
    if ($idx["rf_rate"] != "NA") {
      rf[key]+=$idx["rf_rate"]+0
      rf_count[key]++
    }
    weight[key]+=$idx["weight_calc_ms"]+0
    score[key]+=$idx["optimal_triplet_score"]+0
  }
  END {
    for (key in total) {
      split(key, parts, ",")
      n=total[key]
      avg_gpu = (gpu_count[key] > 0) ? sprintf("%.3f", gpu[key]/gpu_count[key]) : "NA"
      avg_rf = (rf_count[key] > 0) ? sprintf("%.6f", rf[key]/rf_count[key]) : "NA"
      printf "%s,%s,%s,%s,%s,%s,%d,%d,%d,%s,%.3f,%.6f,%.3f,%s,%.3f,%.6f\n",
        parts[1], parts[2], parts[3], parts[4], parts[5], parts[6],
        total[key], success[key]+0, failed[key]+0,
        avg_rf, run[key]/n, stelar[key]/n, cpu[key]/n, avg_gpu, weight[key]/n, score[key]/n
    }
  }
' "$MERGED_CSV" | sort -t, -k4,4n -k3,3 -k5,5 > "$SUMMARY_CSV"

echo "Merged rows: $rows"
echo "Merged CSV:  $MERGED_CSV"
echo "Summary CSV: $SUMMARY_CSV"
