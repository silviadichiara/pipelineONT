#!/bin/bash
set -euo pipefail
source /data/silvia/miniconda3/etc/profile.d/conda.sh
echo "=== USER INPUT PARAMETERS ==="
read -p "Enter the path to the folder containing raw FASTQ files: " REMOTE_DIR
read -p "Enter the path where you want to save the results: " LOCAL_OUTPUT_DIR
read -p "Enter the path where to create the Kraken2 database: " DB_PATH
read -p "Enter the TaxID to filter for (e.g., 10359): " TAXID
read -p "Enter the genome file (after polishing): " GENOME
read -p "Enter the number of threads to use: " THREADS

mkdir -p "$LOCAL_OUTPUT_DIR"
##############################
### 1. KRAKEN2 DB BUILDING ###
##############################
echo "[STEP 1] Checking Kraken2 database..."

if [[ -f "$DB_PATH/hash.k2d" ]]; then
  echo "Kraken2 database already exists at $DB_PATH. Skipping build."
else
  echo "Building Kraken2 database..."
  conda activate kraken2_env
  kraken2-build --standard --db "$DB_PATH" --threads "$THREADS"
  kraken2-build --build --db "$DB_PATH" --threads "$THREADS"
  conda deactivate
fi


####################################
### 2. KRAKEN2 CLASSIFICATION    ###
####################################
echo "[STEP 2] Classifying reads with Kraken2..."
conda activate kraken2_env
mkdir -p "$LOCAL_OUTPUT_DIR/1_kraken"

for R1 in "$REMOTE_DIR"/*_R1.fastq.gz; do
  R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
  SAMPLE=$(basename "$R1" _R1.fastq.gz)

  kraken2 --db "$DB_PATH" --threads "$THREADS" --paired \
    --report "$LOCAL_OUTPUT_DIR/1_kraken/${SAMPLE}.report" \
    --output "$LOCAL_OUTPUT_DIR/1_kraken/${SAMPLE}.kraken" \
    "$R1" "$R2"
done
conda deactivate

######################################
### 3. FILTERING READS WITH SEQTK  ###
######################################
echo "[STEP 3] Filtering reads classified under TaxID $TAXID..."
conda activate kraken2_env  # also used for seqtk
cd "$LOCAL_OUTPUT_DIR/1_kraken"
mkdir -p "$LOCAL_OUTPUT_DIR/2_seqtk_filt"

for KRAKEN_FILE in *.kraken; do
  SAMPLE=$(basename "$KRAKEN_FILE" .kraken)
  R1="$REMOTE_DIR/${SAMPLE}_R1.fastq.gz"
  R2="$REMOTE_DIR/${SAMPLE}_R2.fastq.gz"

  ID_BASE="$LOCAL_OUTPUT_DIR/2_seqtk_filt/${SAMPLE}_tax${TAXID}_IDs"

  awk -v taxid="$TAXID" '$1 == "C" && $3 == taxid { print $2 }' "$KRAKEN_FILE" | sed 's/[ \t].*//' | sed 's/\/[12]$//' | sort | uniq > "${ID_BASE}_all.txt"
  zcat "$R1" | awk 'NR%4==1 {gsub(/^@/, "", $1); split($1,a," "); split(a[1],b,"/"); print b[1]}' | sort | uniq > "${ID_BASE}_r1.txt"
  zcat "$R2" | awk 'NR%4==1 {gsub(/^@/, "", $1); split($1,a," "); split(a[1],b,"/"); print b[1]}' | sort | uniq > "${ID_BASE}_r2.txt"

  comm -12 "${ID_BASE}_all.txt" "${ID_BASE}_r1.txt" | comm -12 - "${ID_BASE}_r2.txt" > "${ID_BASE}_paired.txt"
  awk '{ print "@" $1 }' "${ID_BASE}_paired.txt" > "${ID_BASE}_paired.grep"

  for END in 1 2; do
    zcat "$REMOTE_DIR/${SAMPLE}_R${END}.fastq.gz" | \
    paste - - - - | \
    awk -F'\t' '{ split($1, h, " "); split(h[1], b, "/"); print "@" b[1] "\t" $0 }' | \
    grep -F -f "${ID_BASE}_paired.grep" | \
    cut -f2- | \
    tr '\t' '\n' | \
    gzip > "$LOCAL_OUTPUT_DIR/2_seqtk_filt/${SAMPLE}_tax${TAXID}_R${END}.fastq.gz"
  done
done
conda deactivate

#################################
### 4. QC WITH FASTQC & FASTP ###
#################################
echo "[STEP 4] Running FastQC and fastp..."
conda activate qc_env
mkdir -p "$LOCAL_OUTPUT_DIR/3_fastqc" "$LOCAL_OUTPUT_DIR/4_fastp"

for R1 in "$LOCAL_OUTPUT_DIR/2_seqtk_filt/"*_R1.fastq.gz; do
  SAMPLE=$(basename "$R1" _R1.fastq.gz)
  R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

  fastqc -t 4 -o "$LOCAL_OUTPUT_DIR/3_fastqc" "$R1" "$R2"

  fastp \
    -i "$R1" -I "$R2" \
    -o "$LOCAL_OUTPUT_DIR/4_fastp/${SAMPLE}_R1.fastq.gz" \
    -O "$LOCAL_OUTPUT_DIR/4_fastp/${SAMPLE}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --thread 4 \
    -j "$LOCAL_OUTPUT_DIR/4_fastp/${SAMPLE}_fastp.json" \
    -h "$LOCAL_OUTPUT_DIR/4_fastp/${SAMPLE}_fastp.html"
done
conda deactivate

##################################
### 5. ALIGNMENT + POLISHING   ###
##################################
echo "[STEP 5] Aligning with BWA and polishing with Pilon..."
conda activate illuminareads_env
cd "$LOCAL_OUTPUT_DIR"
mkdir -p 5_bwamem_pol && cd 5_bwamem_pol

if [[ ! -f "${GENOME}.bwt" ]]; then
  bwa index "$GENOME"
fi

BAM_LIST=()
for R1 in "$LOCAL_OUTPUT_DIR/4_fastp/"*_R1.fastq.gz; do
  SAMPLE=$(basename "$R1" _R1.fastq.gz)
  R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
  BAM="${SAMPLE}.bam"
  bwa mem -t "$THREADS" "$GENOME" "$R1" "$R2" | samtools sort -@ "$THREADS" -o "$BAM"
  BAM_LIST+=("$BAM")
done

samtools merge -@ "$THREADS" merged.bam "${BAM_LIST[@]}"
samtools index merged.bam

pilon --genome "$GENOME" \
  --frags merged.bam \
  --output pilon_polished \
  --threads "$THREADS" \
  --changes --vcf --mindepth 5
conda deactivate

###########################
### 6. FINAL QUAST STEP ###
###########################
echo "[STEP 6] Running QUAST..."
mkdir -p "$LOCAL_OUTPUT_DIR/6_quast"
quast.py -t "$THREADS" -o "$LOCAL_OUTPUT_DIR/6_quast" pilon_polished.fasta

########################################
### 7. MAPPING + QUALIMAP ANALYSIS   ###
########################################
echo "[STEP 7] Mapping with minimap2 and running Qualimap..."

MERGED_READS="$LOCAL_OUTPUT_DIR/7_mapping/merged_reads.fastq"
mkdir -p "$LOCAL_OUTPUT_DIR/7_mapping"
cat "$LOCAL_OUTPUT_DIR/4_fastp/"*_R1.fastq.gz "$LOCAL_OUTPUT_DIR/4_fastp/"*_R2.fastq.gz | gunzip -c > "$MERGED_READS"

cd "$LOCAL_OUTPUT_DIR/7_mapping"
conda activate qc_env
minimap2 -ax map-ont "$LOCAL_OUTPUT_DIR/5_bwamem_pol/pilon_polished.fasta" "$MERGED_READS" -t "$THREADS" > mapped.sam

samtools view -bS mapped.sam | samtools sort -o mapped.sorted.bam
samtools index mapped.sorted.bam

qualimap bamqc \
  -bam mapped.sorted.bam \
  -outdir qualimap_report \
  -outformat PDF:HTML
conda deactivate
rm mapped.sam

echo "[DONE] Full RNAseq pipeline completed successfully."
