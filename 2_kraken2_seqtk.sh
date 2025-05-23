#!/bin/bash

# === USER INPUT ===
read -p "Directory containing FASTQ files (from Dorado): " FASTQ_DIR
read -p "Kraken2 executable path: " KRAKEN2_PATH
read -p "Kraken2 database path: " DB_PATH
read -p "Kraken2 output directory: " KRAKEN2_OUT_DIR
read -p "Filtered output directory: " FILTERED_OUT_DIR
read -p "TaxID to filter for (e.g., 10359 for HCMV): " TAXID
THREADS=16

# === PREPARE OUTPUT DIRS ===
mkdir -p "$KRAKEN2_OUT_DIR"
mkdir -p "$FILTERED_OUT_DIR"

# === STEP 1: Run Kraken2 on all FASTQ files ===
echo "[1/2] Running Kraken2 classification..."

for INPUT_FASTQ in "$FASTQ_DIR"/*.fastq; do
    if [ ! -f "$INPUT_FASTQ" ]; then
        echo "No FASTQ files found in $FASTQ_DIR"
        exit 1
    fi

    BASE_NAME=$(basename "$INPUT_FASTQ" .fastq)
    REPORT_FILE="$KRAKEN2_OUT_DIR/${BASE_NAME}.report"
    CLASSIFIED_FILE="$KRAKEN2_OUT_DIR/${BASE_NAME}_classified.out"

    echo "→ Classifying $INPUT_FASTQ..."
    "$KRAKEN2_PATH" \
        --db "$DB_PATH" \
        --output "$CLASSIFIED_FILE" \
        --report "$REPORT_FILE" \
        --threads "$THREADS" \
        "$INPUT_FASTQ"

    if [ $? -ne 0 ]; then
        echo "Error processing $INPUT_FASTQ with Kraken2."
        exit 1
    fi
done

echo "Kraken2 classification completed."

# === STEP 2: Filter reads matching the desired taxid ===
echo "[2/2] Filtering reads matching TaxID $TAXID..."

for CLASSIFIED_FILE in "$KRAKEN2_OUT_DIR"/*_classified.out; do
    BASE_NAME=$(basename "$CLASSIFIED_FILE" _classified.out)
    FASTQ_FILE="$FASTQ_DIR/${BASE_NAME}.fastq"
    ID_LIST="$FILTERED_OUT_DIR/${BASE_NAME}_selected_ids.txt"
    OUTPUT_FASTQ="$FILTERED_OUT_DIR/${BASE_NAME}_taxid${TAXID}.fastq"

    if [ ! -f "$FASTQ_FILE" ]; then
        echo "Corresponding FASTQ not found for $BASE_NAME, skipping."
        continue
    fi

    # Extract read IDs classified as the desired TaxID
    awk -v taxid="$TAXID" '$1 == "C" && $3 == taxid {print $2}' "$CLASSIFIED_FILE" > "$ID_LIST"

    if [ -s "$ID_LIST" ]; then
        seqtk subseq "$FASTQ_FILE" "$ID_LIST" > "$OUTPUT_FASTQ"
        echo "Filtered reads saved to $OUTPUT_FASTQ"
    else
        echo "No reads found for TaxID $TAXID in $BASE_NAME"
    fi
done

echo "All steps completed successfully."
