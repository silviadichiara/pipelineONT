#!/bin/bash
set -euo pipefail

# === Conda setup ===
eval "$(conda shell.bash hook)"
conda activate preprocessing_env

# === USER INPUT ===
read -p "Do you want to convert .fast5 to .pod5? (y/n): " CONVERT_CHOICE

if [[ "$CONVERT_CHOICE" =~ ^[Yy]$ ]]; then
    read -p "Directory containing .fast5 files: " FAST5_DIR
fi

read -p "Directory containing .pod5 files (output or existing): " POD5_DIR
read -p "Output directory for FASTQ files: " OUT_DIR
read -p "Output directory for Dorado basecalling models download: " DORADO_DIR

THREADS=8
LOG_DIR="${OUT_DIR}/logs"

# === SETUP ===
mkdir -p "$POD5_DIR" "$OUT_DIR" "$LOG_DIR" "$DORADO_DIR"

# === [1/3] OPTIONAL: CONVERT FAST5 → POD5 ===
if [[ "$CONVERT_CHOICE" =~ ^[Yy]$ ]]; then
    echo "[1/3] Converting .fast5 to .pod5..."

    shopt -s nullglob
    FAST5_FILES=("$FAST5_DIR"/*.fast5)
    shopt -u nullglob

    if [ ${#FAST5_FILES[@]} -eq 0 ]; then
        echo "No .fast5 files found in $FAST5_DIR. Aborting."
        exit 1
    fi

    for file in "${FAST5_FILES[@]}"; do
        base_name=$(basename "$file" .fast5)
        output_file="${POD5_DIR}/${base_name}.pod5"
        echo "→ Converting $file to $output_file..."
        pod5 convert from_fast5 "$file" -o "$output_file"
    done

    echo "Conversion complete: ${#FAST5_FILES[@]} files processed."
else
    echo "[1/3] Skipping conversion step. Using existing .pod5 files."
fi

# === [2/3] DORADO BASECALLING ===
echo "[2/3] Running Dorado basecalling on each .pod5 file..."

pushd "$DORADO_DIR" > /dev/null
dorado download --model all
echo "Available Dorado basecalling models:"
ls -1
read -p "Name of Dorado basecalling model (e.g. dna_r10.4.1_e8.2_400bps_sup@v4.3.0): " MODEL
popd > /dev/null

shopt -s nullglob
POD5_FILES=("$POD5_DIR"/*.pod5)
shopt -u nullglob

if [ ${#POD5_FILES[@]} -eq 0 ]; then
    echo "No .pod5 files found in $POD5_DIR. Aborting."
    exit 1
fi

for pod5_file in "${POD5_FILES[@]}"; do
    base_name=$(basename "$pod5_file" .pod5)
    output_fastq="${OUT_DIR}/${base_name}.fastq"
    log_file="${LOG_DIR}/basecalling_${base_name}.log"

    echo "→ Basecalling $pod5_file..."
    dorado basecaller "$MODEL" "$pod5_file" \
        --emit-fastq \
        --trim adapters > "$output_fastq" 2> "$log_file"
done

conda deactivate

# === [3/3] DONE ===
echo "[3/3] All tasks completed successfully."
echo "FASTQ output available in: $OUT_DIR"
