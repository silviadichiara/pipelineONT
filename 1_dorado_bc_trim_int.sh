#!/bin/bash

# === USER INPUT ===
read -p "Directory containing .fast5 files: " FAST5_DIR
read -p "Output directory for .pod5 files: " POD5_DIR
read -p "Output directory for FASTQ files: " OUT_DIR
read -p "Path to Dorado basecalling model: " MODEL
THREADS=8

LOG_DIR="$OUT_DIR/logs"

# === SETUP ===
mkdir -p "$POD5_DIR"
mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

# === [1/3] CONVERT FAST5 → POD5 ===
echo "[1/3] Converting .fast5 to .pod5..."

FAST5_COUNT=0
for file in "$FAST5_DIR"/*.fast5; do
    if [ -f "$file" ]; then
        base_name=$(basename "$file" .fast5)
        output_file="$POD5_DIR/${base_name}.pod5"
        
        echo "→ Converting $file to $output_file..."
        pod5 convert from_fast5 "$file" -o "$output_file"

        if [ $? -ne 0 ]; then
            echo "Error converting $file. Aborting."
            exit 1
        fi

        FAST5_COUNT=$((FAST5_COUNT + 1))
    fi
done

if [ "$FAST5_COUNT" -eq 0 ]; then
    echo "No .fast5 files found in $FAST5_DIR. Aborting."
    exit 1
fi

echo "Conversion complete: $FAST5_COUNT files processed."

# === [2/3] DORADO BASECALLING FOR EACH .POD5 FILE ===
echo "[2/3] Running Dorado basecalling on each .pod5 file..."

POD5_COUNT=0
for pod5_file in "$POD5_DIR"/*.pod5; do
    if [ -f "$pod5_file" ]; then
        base_name=$(basename "$pod5_file" .pod5)
        output_fastq="$OUT_DIR/${base_name}.fastq"
        log_file="$LOG_DIR/basecalling_${base_name}.log"

        echo "→ Basecalling $pod5_file..."
        dorado basecaller "$MODEL" "$pod5_file" \
            --emit-fastq \
            --trim adapters > "$output_fastq" 2> "$log_file"

        if [ $? -ne 0 ]; then
            echo "Error during basecalling of $pod5_file. See log: $log_file"
            exit 1
        fi

        POD5_COUNT=$((POD5_COUNT + 1))
    fi
done

if [ "$POD5_COUNT" -eq 0 ]; then
    echo "No .pod5 files found in $POD5_DIR. Aborting."
    exit 1
fi

echo "Basecalling complete for $POD5_COUNT files."

# === [3/3] DONE ===
echo "[3/3] All tasks completed successfully. FASTQ output available in: $OUT_DIR"
