#!/bin/bash
source ~/.bashrc
set -euo pipefail

# === USER INPUT ===
read -p "==> Input reads (FASTQ file): " READS
read -p "==> Output directory: " ROOT
read -p "==> Estimated genome size (e.g. 235k): " GENOME_SIZE
read -p "==> Medaka model (e.g. r1041_e82_400bps_sup_variant_v4.2.0): " MODEL
read -p "==> Threads to use: " THREADS
read -p "==> Project name (e.g. HCMV): " PROJECT_NAME

# === ENVIRONMENT SETUP ===
ASSEMBLY_ENV="assembly_env"
POLISHING_ENV="polishing_env"
PREPROCESSING_ENV="preprocessing_env"
QC_ENV="qc_env"

mkdir -p "$ROOT"

if [[ ! -s "$READS" ]]; then
  echo "[ERROR] Reads file not found or empty: $READS"
  exit 1
fi

echo "[INFO] Activating $ASSEMBLY_ENV..."
conda activate "$ASSEMBLY_ENV"

# === ASSEMBLY TOOLS ===
# NextDenovo
ND_DIR="$ROOT/nextdenovo"
if [[ ! -s "$ND_DIR/nd.asm.fasta" ]]; then
  echo "[STEP] Running NextDenovo..."
  mkdir -p "$ND_DIR"
  echo "$READS" > "$ND_DIR/input.fofn"
  cat > "$ND_DIR/run.cfg" <<EOF
[General]
job_type = local
job_prefix = $PROJECT_NAME
task = all
rewrite = yes
dataroot = $ROOT
workdir = $ND_DIR
input_fofn = input.fofn
read_type = ont

[correct_option]
read_cutoff = 1k
genome_size = $GENOME_SIZE
seed_cutoff = 10k
minimap2_options = -t $THREADS

[assemble_option]
minimap2_options = -t $THREADS
EOF
  nextDenovo "$ND_DIR/run.cfg" | tee "$ND_DIR/nextdenovo.log"
fi

# Canu
Canu_DIR="$ROOT/canu"
Canu_OUT="$Canu_DIR/${PROJECT_NAME}.contigs.fasta"
if [[ ! -s "$Canu_OUT" ]]; then
  echo "[STEP] Running Canu..."
  mkdir -p "$Canu_DIR"
  canu -p "$PROJECT_NAME" -d "$Canu_DIR" genomeSize="$GENOME_SIZE" \
       -nanopore-raw "$READS" useGrid=false minReadLength=500 maxThreads=$THREADS
fi

# Flye
Flye_DIR="$ROOT/flye"
Flye_OUT="$Flye_DIR/assembly.fasta"
if [[ ! -s "$Flye_OUT" ]]; then
  echo "[STEP] Running Flye..."
  mkdir -p "$Flye_DIR"
  flye --nano-corr "$READS" --out-dir "$Flye_DIR" \
       --threads "$THREADS" --genome-size "$GENOME_SIZE" \
       --min-overlap 1000 --asm-coverage 50 | tee "$Flye_DIR/flye.log"
fi

# Wtdbg2
Wtdbg2_DIR="$ROOT/wtdbg2"
Wtdbg2_OUT="$Wtdbg2_DIR/wtdbg_assembly.fasta"
if [[ ! -s "$Wtdbg2_OUT" ]]; then
  echo "[STEP] Running Wtdbg2..."
  mkdir -p "$Wtdbg2_DIR"
  gunzip -c "$Canu_DIR/${PROJECT_NAME}.trimmedReads.fasta.gz" > "$Wtdbg2_DIR/corrected_reads.fasta"
  wtdbg2 -i "$Wtdbg2_DIR/corrected_reads.fasta" -o "$Wtdbg2_DIR/wtdbg_assembly" \
         -t $THREADS -g $GENOME_SIZE -p 0 -k 15 -AS 2 -s 0.05 -L 5000
  wtpoa-cns -t $THREADS -i "$Wtdbg2_DIR/wtdbg_assembly.ctg.lay.gz" -fo "$Wtdbg2_OUT"
fi

# Hifiasm
Hifiasm_DIR="$ROOT/hifiasm"
Hifiasm_OUT="$Hifiasm_DIR/${PROJECT_NAME}_from_gfa.fasta"
if [[ ! -s "$Hifiasm_OUT" ]]; then
  echo "[STEP] Running Hifiasm..."
  mkdir -p "$Hifiasm_DIR"
  hifiasm -o "$Hifiasm_DIR/$PROJECT_NAME" -t $THREADS --ont "$READS" \
      2> "$Hifiasm_DIR/hifiasm_error.log"
  awk '/^S/{print ">"$2"\\n"$3}' "$Hifiasm_DIR/${PROJECT_NAME}.bp.p_ctg.gfa" > "$Hifiasm_OUT"
fi

conda deactivate
conda activate "$POLISHING_ENV"

# === Map of assemblies ===
declare -A assemblies=(
  [canu]="$Canu_OUT"
  [nextdenovo]="$(find "$ND_DIR" -name '*nd.asm*.fasta' | head -n1)"
  [flye]="$Flye_OUT"
  [wtdbg2]="$Wtdbg2_OUT"
  [hifiasm]="$Hifiasm_OUT"
)

# === QUAST before polishing ===
QUAST_PRE="$ROOT/quast_results/pre_polishing"
if [[ ! -d "$QUAST_PRE" ]]; then
  echo "[STEP] Running QUAST (pre-polishing)..."
  mkdir -p "$QUAST_PRE"
  quast_cmd="quast.py -t $THREADS -o $QUAST_PRE"
  for tool in "${!assemblies[@]}"; do
    [[ -s "${assemblies[$tool]}" ]] && quast_cmd+=" ${assemblies[$tool]}"
  done
  eval "$quast_cmd"
fi

# === POLISHING LOOP (Racon + Medaka) ===
POLISH_DIR="$ROOT/polished_assemblies"
mkdir -p "$POLISH_DIR"
for tool in "${!assemblies[@]}"; do
  ASSEMBLY="${assemblies[$tool]}"
  [[ ! -s "$ASSEMBLY" ]] && continue

  OUT_DIR="$POLISH_DIR/${tool}_polish"
  if [[ -s "$OUT_DIR/medaka/consensus.fasta" ]]; then
    echo "[INFO] Skipping polishing for $tool (already done)."
    continue
  fi

  echo "[STEP] Polishing $tool..."
  mkdir -p "$OUT_DIR"

  conda activate "$PREPROCESSING_ENV"
  minimap2 -ax map-ont -t $THREADS "$ASSEMBLY" "$READS" > "$OUT_DIR/reads.sam"
  conda deactivate

  conda activate "$POLISHING_ENV"
  racon -t $THREADS "$READS" "$OUT_DIR/reads.sam" "$ASSEMBLY" > "$OUT_DIR/racon_polished.fasta"
  medaka_consensus -i "$READS" -d "$OUT_DIR/racon_polished.fasta" -o "$OUT_DIR/medaka" -m "$MODEL"
done

# === QUAST after polishing ===
QUAST_POST="$ROOT/quast_results/post_polishing"
if [[ ! -d "$QUAST_POST" ]]; then
  echo "[STEP] Running QUAST (post-polishing)..."
  mkdir -p "$QUAST_POST"
  POST_FILES=$(find "$POLISH_DIR" -name "consensus.fasta")
  quast.py -t $THREADS -o "$QUAST_POST" $POST_FILES
fi

conda deactivate

# === MAPPING + QUALIMAP ===
MAPPING_DIR="$ROOT/mapping"
QUALIMAP_DIR="$MAPPING_DIR/qualimap_reports"
mkdir -p "$MAPPING_DIR" "$QUALIMAP_DIR"

for fasta in "$POLISH_DIR"/*/medaka/consensus.fasta; do
  tool=$(basename "$(dirname "$(dirname "$fasta")")")
  sam="$MAPPING_DIR/${tool}.sam"
  bam="$MAPPING_DIR/${tool}.sorted.bam"

  echo "[STEP] Mapping reads on $tool..."
  conda activate "$PREPROCESSING_ENV"
  minimap2 -ax map-ont "$fasta" "$READS" > "$sam"
  conda deactivate

  conda activate "$QC_ENV"
  samtools view -bS "$sam" | samtools sort -o "$bam"
  samtools index "$bam"
  qualimap bamqc -bam "$bam" -outdir "$QUALIMAP_DIR/${tool}_qualimap" -outformat PDF:HTML --java-mem-size=4G
  rm "$sam"
  conda deactivate
done

echo "Pipeline completed. Reports and polished assemblies are in: $ROOT"
echo "Before continuing, have a look to the denovo assemblies' stats and choose the best assembler (e.g. Flye)"
