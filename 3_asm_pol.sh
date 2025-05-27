#!/bin/bash
source ~/.bashrc
set -euo pipefail

# === USER INPUT ===
read -p "==> Input reads (FASTQ file): " READS
read -p "==> Output directory: " ROOT
read -p "==> Estimated genome size (e.g. 235k): " GENOME_SIZE
read -p "==> Medaka model (e.g. r1041_e82_400bps_sup_variant_v4.2.0): " MODEL

ASSEMBLY_ENV="assembly_env"
POLISHING_ENV="polishing_env"
THREADS=8
PROJECT_NAME="HCMV"

mkdir -p "$ROOT"

if [ ! -s "$READS" ]; then
  echo "Error: reads file not found or empty: $READS"
  exit 1
fi

conda activate "$ASSEMBLY_ENV"

# NEXTDENOVO
ND_DIR="$ROOT/nextdenovo"
if [ -d "$ND_DIR" ] && [ "$(ls -A "$ND_DIR")" ]; then
  echo "Skipping NextDenovo: already completed."
else
  echo "Running NextDenovo..."
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

# CANU
Canu_DIR="$ROOT/canu"
if [ -d "$Canu_DIR" ] && [ "$(ls -A "$Canu_DIR")" ]; then
  echo "Skipping Canu: already completed."
else
  echo "Running Canu..."
  mkdir -p "$Canu_DIR"
  canu -p "$PROJECT_NAME" -d "$Canu_DIR" genomeSize="$GENOME_SIZE" -nanopore-raw "$READS" useGrid=false minReadLength=500 maxThreads=$THREADS
fi

# FLYE
Flye_DIR="$ROOT/flye"
if [ -d "$Flye_DIR" ] && [ "$(ls -A "$Flye_DIR")" ]; then
  echo "Skipping Flye: already completed."
else
  echo "Running Flye..."
  mkdir -p "$Flye_DIR"
  flye --nano-corr "$READS" --out-dir "$Flye_DIR" --threads "$THREADS" --genome-size "$GENOME_SIZE" --min-overlap 1000 --asm-coverage 50 | tee "$Flye_DIR/flye.log"
fi

# WTDBG2
Wtdbg2_DIR="$ROOT/wtdbg2"
if [ -d "$Wtdbg2_DIR" ] && [ "$(ls -A "$Wtdbg2_DIR")" ]; then
  echo "Skipping WTDBG2: already completed."
else
  echo "Running Wtdbg2..."
  mkdir -p "$Wtdbg2_DIR"
  gunzip -c "$Canu_DIR/${PROJECT_NAME}.trimmedReads.fasta.gz" > "$Wtdbg2_DIR/corrected_reads.fasta"
  wtdbg2 -i "$Wtdbg2_DIR/corrected_reads.fasta" -o "$Wtdbg2_DIR/wtdbg_assembly" -t $THREADS -g $GENOME_SIZE -p 0 -k 15 -AS 2 -s 0.05 -L 5000
  wtpoa-cns -t $THREADS -i "$Wtdbg2_DIR/wtdbg_assembly.ctg.lay.gz" -fo "$Wtdbg2_DIR/wtdbg_assembly.fasta"
fi

# HIFIASM
Hifiasm_DIR="$ROOT/hifiasm"
HIFIASM_FASTA="$Hifiasm_DIR/${PROJECT_NAME}_from_gfa.fasta"
if [ -s "$HIFIASM_FASTA" ]; then
  echo "Skipping Hifiasm: already completed."
else
  echo "Running Hifiasm..."
  mkdir -p "$Hifiasm_DIR"
  hifiasm -o "$Hifiasm_DIR/$PROJECT_NAME" -t $THREADS --ont "$READS" 2> "$Hifiasm_DIR/hifiasm_error.log"
  awk '/^S/{print ">"$2"\n"$3}' "$Hifiasm_DIR/${PROJECT_NAME}.bp.p_ctg.gfa" > "$HIFIASM_FASTA"
fi

conda deactivate
conda activate "$POLISHING_ENV"

 declare -A assemblies=(
    ["canu"]="$Canu_DIR/${PROJECT_NAME}.contigs.fasta"
    ["nextdenovo"]="$(find "$ND_DIR" -name '*nd.asm*.fasta' | head -n1)"
    ["flye"]="$Flye_DIR/assembly.fasta"
    ["wtdbg2"]="$Wtdbg2_DIR/wtdbg_assembly.fasta"
    ["hifiasm"]="$HIFIASM_FASTA"
  )
# QUAST PRE
QUAST_PRE="$ROOT/quast_results/pre_polishing"
if [ -d "$QUAST_PRE" ] && [ "$(ls -A "$QUAST_PRE")" ]; then
  echo "Skipping QUAST pre-polishing: already completed."
else
  echo "Running QUAST pre-polishing..."
  mkdir -p "$QUAST_PRE"
  quast_cmd="quast.py -t $THREADS -o $QUAST_PRE"
  for tool in "${!assemblies[@]}"; do
    [ -s "${assemblies[$tool]}" ] && quast_cmd+=" ${assemblies[$tool]}"
  done
  eval "$quast_cmd"
fi

# POLISHING: Racon + Medaka
POLISH_DIR="$ROOT/polished_assemblies"
mkdir -p "$POLISH_DIR"
for ASSEMBLY in "${assemblies[@]}"; do
  [ ! -s "$ASSEMBLY" ] && continue
  SAMPLE=$(basename "$ASSEMBLY" .fasta)
  OUT_DIR="$POLISH_DIR/$SAMPLE"
  if [ -s "$OUT_DIR/medaka/consensus.fasta" ]; then
    echo "Skipping polishing for $SAMPLE: already completed."
    continue
  fi
  echo "--> Polishing $SAMPLE"
  mkdir -p "$OUT_DIR"
  minimap2 -ax map-ont -t $THREADS "$ASSEMBLY" "$READS" > "$OUT_DIR/reads.sam"
  racon -t $THREADS "$READS" "$OUT_DIR/reads.sam" "$ASSEMBLY" > "$OUT_DIR/racon_polished.fasta"
  medaka_consensus -i "$READS" -d "$OUT_DIR/racon_polished.fasta" -o "$OUT_DIR/medaka" -m "$MODEL"

done

# QUAST POST
QUAST_POST="$ROOT/quast_results/post_polishing"
if [ -d "$QUAST_POST" ] && [ "$(ls -A "$QUAST_POST")" ]; then
  echo "Skipping QUAST post-polishing: already completed."
else
  echo "Running QUAST post-polishing..."
  mkdir -p "$QUAST_POST"
  post_files=$(find "$POLISH_DIR" -name "consensus.fasta")
  quast.py -t $THREADS -o "$QUAST_POST" $post_files
fi

conda deactivate

echo "Pipeline completed."
