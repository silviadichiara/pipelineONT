#!/bin/bash

# Check for glsearch36
if ! command -v glsearch36 &> /dev/null; then
  echo "glsearch36 command not found. Please install it."
  exit 1
fi

# Prompt user input
read -p "Enter the path to the input masked FASTA file: " input_file
read -p "Enter the path to the output unaligned regions FASTA file (.fasta): " output_file
read -p "Enter the path to the output alignment directory: " alignment_dir
read -p "Enter the path to the reference database FASTA file: " db_file

# Check input files
if [[ ! -f "$input_file" ]]; then
  echo "Masked input FASTA file not found!"
  exit 1
fi

if [[ ! -f "$db_file" ]]; then
  echo "Reference database file not found!"
  exit 1
fi

mkdir -p "$alignment_dir"
report_file="${alignment_dir}/alignment_report.tsv"
echo -e "Region_ID\tBest_Hit_ID\t%_Identity\tScore\tQuery_Pos\tCoverage(%)" > "$report_file"

header=""
sequence=""

while read -r line; do
  if [[ "$line" =~ ^\> ]]; then
    header="${line#>}"
  else
    sequence+="$line"
  fi
done < "$input_file"

seq_len=${#sequence}
start=0
in_region=0
region_count=1
> "$output_file"

process_region() {
  local s=$1
  local e=$2
  local subseq="${sequence:$((s-1)):$((e - s + 1))}"
  local region_id="${header}_region_${region_count}_${s}_${e}"

  echo ">${region_id}" >> "$output_file"
  echo "$subseq" >> "$output_file"

  region_file="${alignment_dir}/${region_id}.fasta"
  echo ">${region_id}" > "$region_file"
  echo "$subseq" >> "$region_file"

  gl_out="${alignment_dir}/${region_id}_glsearch.out"
  glsearch36 -m 8c "$region_file" "$db_file" > "$gl_out"

  hit_line=$(grep -v '^#' "$gl_out" | head -n 1)
  if [[ -n "$hit_line" ]]; then
    read -r _ best_hit_id percent_identity _ _ _ q_start q_end _ _ _ score _ <<< "$hit_line"
    query_len=${#subseq}
    aligned_len=$(( q_end > q_start ? q_end - q_start + 1 : q_start - q_end + 1 ))
    coverage=$(awk "BEGIN { printf \"%.2f\", ($aligned_len / $query_len) * 100 }")
    echo -e "${region_id}\t${best_hit_id}\t${percent_identity}\t${score}\t${q_start}-${q_end}\t${coverage}" >> "$report_file"
  else
    echo -e "${region_id}\tNO_HIT\t-\t-\t-\t-" >> "$report_file"
  fi

  region_count=$((region_count + 1))
}

for (( i=0; i<seq_len; i++ )); do
  char="${sequence:$i:1}"
  pos=$((i+1))
  if [[ "$char" != "N" && $in_region -eq 0 ]]; then
    start=$pos
    in_region=1
  elif ([[ "$char" == "N" ]] || [[ $pos -eq $seq_len ]]) && [[ $in_region -eq 1 ]]; then
    end=$(( pos == seq_len && char != "N" ? pos : pos - 1 ))
    process_region "$start" "$end"
    in_region=0
  fi
done

# Handle last region if sequence ends without N
if [[ $in_region -eq 1 ]]; then
  process_region "$start" "$seq_len"
fi

echo "All unmasked regions extracted to: $output_file"
echo "Alignments stored in: $alignment_dir"
echo "Alignment summary report saved to: $report_file"
