#!/bin/bash
set -euo pipefail
shopt -s nullglob

# === LOAD CONDA (portable) ===
if ! command -v conda &> /dev/null; then
  echo "[ERROR] Conda not found. Please ensure conda is installed and initialized."
  exit 1
fi
eval "$(conda shell.bash hook)"

# === FUNCTION: Check if step should run ===
should_run_step() {
  local output_check="$1"
  local step_name="$2"
  local step_dir
  step_dir=$(dirname "$output_check")

  if [[ -e "$output_check" ]]; then
    echo "[INFO] Output for '$step_name' already exists at: $output_check"
    read -p "Do you want to skip this step? (y to skip / n to overwrite): " choice
    if [[ "$choice" =~ ^[Yy]$ ]]; then
      echo "[INFO] Skipping $step_name..."
      return 1
    else
      echo "[INFO] Overwriting $step_name: removing $step_dir..."
      rm -rf "$step_dir"
      mkdir -p "$step_dir"
      return 0
    fi
  else
    return 0
  fi
}

# === USER INPUT ===
read -p "Enter ONT reads FASTQ file path: " READS
read -p "Enter reference genome FASTA file path: " REFERENCE
read -p "Enter de novo genome FASTA file path: " DENOVO_GENOME
read -p "Enter annotation GFF3 file path (dowload from NCBI): " GFF
read -p "Enter output directory path: " OUTDIR
read -p "Enter number of threads to use: " THREADS
read -p "Enter Medaka model to use: " MODEL

# === CHECK FILES ===
for FILE in "$READS" "$REFERENCE" "$DENOVO_GENOME" "$GFF"; do
  [[ -f "$FILE" ]] || { echo "[ERROR] Missing input: $FILE"; exit 1; }
done
mkdir -p "$OUTDIR"

# === STEP 0: MUMmer Comparison ===
if should_run_step "$OUTDIR/mummer/ref_vs_denovo.coords" "MUMmer Comparison"; then
  echo "[0] Running MUMmer comparison..."
  mkdir -p "$OUTDIR/mummer"
  conda activate mummer_env
  nucmer --prefix="$OUTDIR/mummer/ref_vs_denovo" "$REFERENCE" "$DENOVO_GENOME"
  delta-filter -1 "$OUTDIR/mummer/ref_vs_denovo.delta" > "$OUTDIR/mummer/ref_vs_denovo.filtered.delta"
  show-coords -rcl "$OUTDIR/mummer/ref_vs_denovo.filtered.delta" > "$OUTDIR/mummer/ref_vs_denovo.coords"
  conda deactivate
fi

# === STEP 1: De novo vs Reference Variant Calling ===
if should_run_step "$OUTDIR/denovo_vs_ref/denovo_vs_ref.vcf.gz" "De novo Variant Calling"; then
  echo "[1] Variant calling from de novo assembly..."
  mkdir -p "$OUTDIR/denovo_vs_ref"
  conda activate preprocessing_env
  minimap2 -ax asm5 "$REFERENCE" "$DENOVO_GENOME" > "$OUTDIR/denovo_vs_ref/aln.sam"
  conda deactivate

  conda activate qc_env
  samtools sort -o "$OUTDIR/denovo_vs_ref/aln.sorted.bam" "$OUTDIR/denovo_vs_ref/aln.sam"
  samtools index "$OUTDIR/denovo_vs_ref/aln.sorted.bam"
  rm "$OUTDIR/denovo_vs_ref/aln.sam"
  bcftools mpileup -f "$REFERENCE" "$OUTDIR/denovo_vs_ref/aln.sorted.bam" | \
    bcftools call -mv -Oz -o "$OUTDIR/denovo_vs_ref/denovo_vs_ref.vcf.gz"
  bcftools index "$OUTDIR/denovo_vs_ref/denovo_vs_ref.vcf.gz"
  conda deactivate
fi

# === STEP 2: Medaka Consensus + Variant Calling ===
if should_run_step "$OUTDIR/medaka/variant_calling.vcf.gz" "Medaka Variant Calling"; then
  echo "[2] Medaka consensus and variant calling..."
  mkdir -p "$OUTDIR/medaka"

  conda activate preprocessing_env
  minimap2 -ax map-ont "$REFERENCE" "$READS" > "$OUTDIR/medaka/aln.sam"
  conda deactivate

  conda activate qc_env
  samtools sort -o "$OUTDIR/medaka/aln.sorted.bam" "$OUTDIR/medaka/aln.sam"
  samtools index "$OUTDIR/medaka/aln.sorted.bam"
  rm "$OUTDIR/medaka/aln.sam"
  conda deactivate

  conda activate polishing_env
  CONS_HDF="$OUTDIR/medaka/medaka_consensus.hdf"
  medaka consensus "$OUTDIR/medaka/aln.sorted.bam" "$CONS_HDF" --threads "$THREADS" --model "$MODEL"
  [[ -s "$CONS_HDF" ]] || { echo "[ERROR] Medaka output missing."; exit 1; }

  medaka variant "$REFERENCE" "$CONS_HDF" "$OUTDIR/medaka/variant_calling.vcf"
  bgzip -c "$OUTDIR/medaka/variant_calling.vcf" > "$OUTDIR/medaka/variant_calling.vcf.gz"
  bcftools index "$OUTDIR/medaka/variant_calling.vcf.gz"
  conda deactivate
fi

# === STEP 3: Compare VCFs ===
if should_run_step "$OUTDIR/final/vcf_comparison/sites.txt" "VCF Comparison"; then
  echo "[3] Comparing Medaka and de novo VCFs..."
  mkdir -p "$OUTDIR/final/vcf_comparison"
  conda activate qc_env
  bcftools isec -p "$OUTDIR/final/vcf_comparison" \
    "$OUTDIR/medaka/variant_calling.vcf.gz" \
    "$OUTDIR/denovo_vs_ref/denovo_vs_ref.vcf.gz"
  conda deactivate
fi

# === STEP 4: Filter Medaka VCF ===
if should_run_step "$OUTDIR/final/highconf.filtered.vcf.gz" "Medaka VCF Filtering"; then
  echo "[4] Filtering Medaka VCF (QUAL >= 30)..."
  VCF_FILTERED="$OUTDIR/final/highconf.filtered.vcf"
  VCF_COMPRESSED="$VCF_FILTERED.gz"
  bcftools filter -e 'QUAL<30' "$OUTDIR/medaka/variant_calling.vcf.gz" -o "$VCF_FILTERED" -Ov
  bgzip -c "$VCF_FILTERED" > "$VCF_COMPRESSED"
  bcftools index "$VCF_COMPRESSED"
fi

# === STEP 5: Generate consensus from filtered VCF ===
if should_run_step "$OUTDIR/final/consensus_filtered.fasta" "Consensus Generation"; then
  echo "[5] Generating filtered consensus genome..."
  bcftools consensus -f "$REFERENCE" "$OUTDIR/final/highconf.filtered.vcf.gz" > "$OUTDIR/final/consensus_filtered.fasta"
fi

echo "[DONE] Variant analysis pipeline completed successfully."
echo "→ Consensus genome: $OUTDIR/final/consensus_filtered.fasta"
