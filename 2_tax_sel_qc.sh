#!/bin/bash
set -euo pipefail

# === USER INPUT ===
read -p "Directory containing FASTQ files (from Dorado): " FASTQ_DIR
read -p "Kraken2 database path: " DB_PATH
read -p "Kraken2 output directory: " KRAKEN2_OUT_DIR
read -p "Filtered output directory: " FILTERED_OUT_DIR
read -p "TaxID to filter for (e.g., 10359 for HCMV): " TAXID
read -p "NanoPlot output directory: " NANOPLOT_DIR

THREADS=8
MERGED_FASTQ="$FILTERED_OUT_DIR/merged_taxid${TAXID}.fastq"

# === SETUP ===
echo "[INFO] Creating output directories..."
mkdir -p "$KRAKEN2_OUT_DIR" "$FILTERED_OUT_DIR" "$NANOPLOT_DIR"

# === [1/4] Kraken2 Classification ===
echo "[1/4] Running Kraken2 classification..."
conda activate kraken2_env
shopt -s nullglob
FASTQ_FILES=("$FASTQ_DIR"/*.fastq)
if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "[ERROR] No FASTQ files found in $FASTQ_DIR"
    exit 1
fi

for INPUT_FILE in "${FASTQ_FILES[@]}"; do
    BASE=$(basename "$INPUT_FILE" .fastq)
    REPORT="$KRAKEN2_OUT_DIR/${BASE}.report"
    CLASSIFIED="$KRAKEN2_OUT_DIR/${BASE}_classified.out"

    echo "→ Classifying $INPUT_FILE..."
    kraken2 --db "$DB_PATH" \
        --output "$CLASSIFIED" \
        --report "$REPORT" \
        --threads "$THREADS" \
        "$INPUT_FILE" || {
            echo "[ERROR] Kraken2 failed on $INPUT_FILE"
            exit 1
        }
done

echo "[INFO] Kraken2 classification complete."

# === [2/4] TaxID Filtering ===
echo "[2/4] Filtering reads for TaxID $TAXID..."
> "$MERGED_FASTQ"  # Reset merge file

for CLASSIFIED_FILE in "$KRAKEN2_OUT_DIR"/*_classified.out; do
    BASE=$(basename "$CLASSIFIED_FILE" _classified.out)
    FASTQ="$FASTQ_DIR/${BASE}.fastq"
    IDS="$FILTERED_OUT_DIR/${BASE}_selected_ids.txt"
    OUT_FASTQ="$FILTERED_OUT_DIR/${BASE}_taxid${TAXID}.fastq"

    if [ ! -f "$FASTQ" ]; then
        echo "[WARNING] No matching FASTQ for $BASE. Skipping."
        continue
    fi

    awk -v taxid="$TAXID" '$1 == "C" && $3 == taxid {print $2}' "$CLASSIFIED_FILE" > "$IDS"

    if [ -s "$IDS" ]; then
        seqtk subseq "$FASTQ" "$IDS" > "$OUT_FASTQ"
        cat "$OUT_FASTQ" >> "$MERGED_FASTQ"
        echo "[INFO] Filtered and added: $OUT_FASTQ"
    else
        echo "[INFO] No reads for TaxID $TAXID in $BASE"
    fi
done

if [ ! -s "$MERGED_FASTQ" ]; then
    echo "[ERROR] No reads retained after filtering. Exiting."
    exit 1
fi
conda deactivate
# === [3/4] NanoPlot QC ===
echo "[3/4] Running NanoPlot QC..."
conda activate qc_env
if ! command -v NanoPlot &> /dev/null; then
    echo "[ERROR] NanoPlot is not installed. Aborting."
    exit 1
fi

NanoPlot \
    --fastq "$MERGED_FASTQ" \
    --loglength \
    --threads "$THREADS" \
    --plots dot kde \
    --title "QC Report - merged_taxid${TAXID}.fastq" \
    --prefix "HCMV_taxid${TAXID}_report" \
    --outdir "$NANOPLOT_DIR"
conda deactivate
# === [4/4] DONE ===
echo "[4/4] NanoPlot complete. Results saved in: $NANOPLOT_DIR"
