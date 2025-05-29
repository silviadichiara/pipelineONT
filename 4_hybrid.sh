#!/bin/bash
set -euo pipefail
shopt -s nullglob
source /data/silvia/miniconda3/etc/profile.d/conda.sh

should_run_step() {
  local output_check="$1"
  local step_name="$2"

  if [[ -e "$output_check" ]]; then
    echo "[INFO] Output for '$step_name' already exists at: $output_check"
    read -p "Do you want to skip this step? (y to skip / n to overwrite): " choice
    if [[ "$choice" =~ ^[Yy]$ ]]; then
      echo "Skipping $step_name..."
      return 1
    else
      echo "Overwriting $step_name..."
      return 0
    fi
  else
    return 0
  fi
}

# === USER INPUT PARAMETERS ===
read -p "Enter the path to the folder containing Illumina reads raw FASTQ files: " REMOTE_DIR
read -p "Enter the path where you want to save the results: " LOCAL_OUTPUT_DIR
read -p "Enter the path where to create the Kraken2 database: " DB_PATH
read -p "Enter the TaxID to filter for (e.g., 10359): " TAXID
read -p "Enter the path to the base folder containing polished assemblies: " BASE_GENOME_PATH
read -p "Enter the assembler used (e.g., flye, canu, etc.): " ASSEMBLER
read -p "Enter the number of threads to use: " THREADS

GENOME="${BASE_GENOME_PATH}/${ASSEMBLER}_polish/medaka/consensus.fasta"
if [[ ! -f "$GENOME" ]]; then
  echo "[ERROR] Polished genome not found at $GENOME"
  exit 1
fi

mkdir -p "$LOCAL_OUTPUT_DIR"
echo "Using polished genome: $GENOME"

# STEP 1: KRAKEN2 DB
if should_run_step "$DB_PATH/hash.k2d" "Kraken2 DB Build"; then
  echo "[STEP 1] Building Kraken2 database..."
  conda activate kraken2_env
  kraken2-build --standard --db "$DB_PATH" --threads "$THREADS"
  kraken2-build --build --db "$DB_PATH" --threads "$THREADS"
  conda deactivate
fi

# STEP 2: KRAKEN2 CLASSIFICATION
if should_run_step "$LOCAL_OUTPUT_DIR/1_kraken" "Kraken2 Classification"; then
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
fi

# STEP 3: FILTERING WITH SEQTK
if should_run_step "$LOCAL_OUTPUT_DIR/2_seqtk_filt" "Read Filtering with seqtk"; then
  echo "[STEP 3] Filtering reads classified under TaxID $TAXID..."
  conda activate kraken2_env
  cd "$LOCAL_OUTPUT_DIR/1_kraken"
  mkdir -p "$LOCAL_OUTPUT_DIR/2_seqtk_filt"
  for KRAKEN_FILE in *.kraken; do
    SAMPLE=$(basename "$KRAKEN_FILE" .kraken)
    R1="$REMOTE_DIR/${SAMPLE}_R1.fastq.gz"
    R2="$REMOTE_DIR/${SAMPLE}_R2.fastq.gz"
    ID_BASE="$LOCAL_OUTPUT_DIR/2_seqtk_filt/${SAMPLE}_tax${TAXID}_IDs"
    awk -v taxid="$TAXID" '$1 == "C" && $3 == taxid { print $2 }' "$KRAKEN_FILE" | sed 's/[ \t].*//' | sed 's\/\/[12]$//' | sort | uniq > "${ID_BASE}_all.txt"
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
fi

# STEP 4: FASTQC & FASTP
if should_run_step "$LOCAL_OUTPUT_DIR/3_fastqc" "FastQC and fastp"; then
  echo "[STEP 4] Running FastQC and fastp..."
  conda activate qc_env
  mkdir -p "$LOCAL_OUTPUT_DIR/3_fastqc" "$LOCAL_OUTPUT_DIR/4_fastp"
  for R1 in "$LOCAL_OUTPUT_DIR"/2_seqtk_filt/*_R1.fastq.gz; do
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
fi

# STEP 5: BWA + PILON
if should_run_step "$LOCAL_OUTPUT_DIR/5_bwamem_pol/pilon_polished.fasta" "BWA and Pilon Polishing"; then
  echo "[STEP 5] Aligning with BWA and polishing with Pilon..."
  conda activate illuminareads_env
  mkdir -p "$LOCAL_OUTPUT_DIR/5_bwamem_pol"
  cd "$LOCAL_OUTPUT_DIR/5_bwamem_pol"

  if [[ ! -f "${GENOME}.bwt" ]]; then
    bwa index "$GENOME"
  fi

  BAM_LIST=()
  for R1 in "$LOCAL_OUTPUT_DIR"/4_fastp/*_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    BAM="${SAMPLE}.bam"
    bwa mem -t "$THREADS" "$GENOME" "$R1" "$R2" | samtools sort -@ "$THREADS" -o "$BAM"
    samtools index "$BAM"
    BAM_LIST+=("$BAM")
  done

  samtools merge -@ "$THREADS" merged.bam "${BAM_LIST[@]}"
  samtools index merged.bam

  read -p "Insert the path to pilon.jar: " PILON_JAR
  if [[ ! -f "$PILON_JAR" ]]; then
    echo "[ERROR] pilon.jar not found in $PILON_JAR"
    exit 1
  fi

  read -p "Insert n. of GB to be assigned to Java (eg. 32G): " JAVA_RAM

  echo "Executing Pilon with $JAVA_RAM RAM using $PILON_JAR..."
  java -Xmx${JAVA_RAM} -jar "$PILON_JAR" \
    --genome "$GENOME" \
    --frags merged.bam \
    --output pilon_polished \
    --changes --vcf --mindepth 5 --fix all \
    > pilon.log 2>&1

  conda deactivate
fi

# STEP 6: QUAST
if should_run_step "$LOCAL_OUTPUT_DIR/6_quast" "QUAST Analysis"; then
  echo "[STEP 6] Running QUAST..."
  mkdir -p "$LOCAL_OUTPUT_DIR/6_quast"
  quast.py -t "$THREADS" -o "$LOCAL_OUTPUT_DIR/6_quast" "$LOCAL_OUTPUT_DIR/5_bwamem_pol/pilon_polished.fasta"
fi

# STEP 7: MAPPING + QUALIMAP
if should_run_step "$LOCAL_OUTPUT_DIR/7_mapping/qualimap_report" "Qualimap Analysis"; then
  echo "[STEP 7] Mapping with minimap2 and running Qualimap..."
  mkdir -p "$LOCAL_OUTPUT_DIR/7_mapping"
  read -p "Enter path to filtered FASTQ files (use wildcard, e.g. /path/*.fastq): " READ_GLOB
  READ_FILES=( $READ_GLOB )
  if [[ ${#READ_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No FASTQ files found matching pattern: $READ_GLOB" >&2
    exit 1
  fi
  MERGED_READS="$LOCAL_OUTPUT_DIR/7_mapping/merged_reads.fastq"
  echo "[INFO] Merging ${#READ_FILES[@]} FASTQ files..."
  cat "${READ_FILES[@]}" > "$MERGED_READS"
  cd "$LOCAL_OUTPUT_DIR/7_mapping"
  conda activate preprocessing_env
  minimap2 -ax map-ont "$LOCAL_OUTPUT_DIR/5_bwamem_pol/pilon_polished.fasta" "$MERGED_READS" -t "$THREADS" > mapped.sam
  conda deactivate
  conda activate qc_env
  samtools view -bS mapped.sam | samtools sort -o mapped.sorted.bam
  samtools index mapped.sorted.bam
  qualimap bamqc \
    -bam mapped.sorted.bam \
    -outdir qualimap_report \
    -outformat PDF:HTML \
    --java-mem-size=16G
  conda deactivate
  rm mapped.sam
fi

echo "[DONE] Full RNAseq pipeline completed successfully."
