#!/bin/bash

# === USER INPUT ===
read -p "Directory containing FASTQ files (from Dorado): " FASTQ_DIR
read -p "Kraken2 executable path: " KRAKEN2_PATH
read -p "Kraken2 database path: " DB_PATH
read -p "Kraken2 output directory: " KRAKEN2_OUT_DIR
read -p "Filtered output directory: " FILTERED_OUT_DIR
read -p "TaxID to filter for (e.g., 10359 for HCMV): " TAXID
read -p "Merged FASTQ file for NanoPlot: " INPUT_FASTQ
read -p "NanoPlot output directory: " NANOPLOT_DIR

THREADS=8

# === PREPARE OUTPUT DIRS ===
mkdir -p "$KRAKEN2_OUT_DIR"
mkdir -p "$FILTERED_OUT_DIR"
mkdir -p "$NANOPLOT_DIR"

# === STEP 1: Kraken2 Classification ===
echo "[1/3] Running Kraken2 classification..."

for INPUT_FILE in "$FASTQ_DIR"/*.fastq; do
    [ -f "$INPUT_FILE" ] || { echo "No FASTQ files found in $FASTQ_DIR"; exit 1; }

    BASE=$(basename "$INPUT_FILE" .fastq)
    REPORT="$KRAKEN2_OUT_DIR/${BASE}.report"
    CLASSIFIED="$KRAKEN2_OUT_DIR/${BASE}_classified.out"

    echo "→ Classifying $INPUT_FILE..."
    "$KRAKEN2_PATH" --db "$DB_PATH" \
        --output "$CLASSIFIED" \
        --report "$REPORT" \
        --threads "$THREADS" \
        "$INPUT_FILE"

    [ $? -eq 0 ] || { echo "Kraken2 failed on $INPUT_FILE"; exit 1; }
done

echo "Kraken2 classification complete."

# === STEP 2: TaxID Filtering ===
echo "[2/3] Filtering reads for TaxID $TAXID..."

for CLASSIFIED_FILE in "$KRAKEN2_OUT_DIR"/*_classified.out; do
    BASE=$(basename "$CLASSIFIED_FILE" _classified.out)
    FASTQ="$FASTQ_DIR/${BASE}.fastq"
    IDS="$FILTERED_OUT_DIR/${BASE}_selected_ids.txt"
    OUT_FASTQ="$FILTERED_OUT_DIR/${BASE}_taxid${TAXID}.fastq"

    [ -f "$FASTQ" ] || { echo "No matching FASTQ for $BASE"; continue; }

    awk -v taxid="$TAXID" '$1 == "C" && $3 == taxid {print $2}' "$CLASSIFIED_FILE" > "$IDS"

    if [ -s "$IDS" ]; then
        seqtk subseq "$FASTQ" "$IDS" > "$OUT_FASTQ"
        echo "Filtered FASTQ: $OUT_FASTQ"
    else
        echo "No matching reads for TaxID $TAXID in $BASE"
    fi
done

# === STEP 3: NanoPlot QC ===
echo "[3/3] Running NanoPlot..."

if ! command -v NanoPlot &> /dev/null; then
    echo "Error: NanoPlot is not installed."
    exit 1
fi

NanoPlot \
    --fastq "$INPUT_FASTQ" \
    --loglength \
    --threads "$THREADS" \
    --plots dot kde \
    --title "QC Report - $(basename "$INPUT_FASTQ")" \
    --prefix HCMV_report \
    --outdir "$NANOPLOT_DIR"

echo "NanoPlot analysis complete. Results in: $NANOPLOT_DIR"
