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

# === CHECK INPUT FILE ===
if [ ! -s "$READS" ]; then
  echo "Error: reads file not found or empty: $READS"
  exit 1
fi

echo "Activating assembly environment: $ASSEMBLY_ENV"
conda init
conda activate "$ASSEMBLY_ENV"

########################################
# 1. NEXTDENOVO
########################################
echo "Running NextDenovo..."
ND_DIR="$ROOT/nextdenovo"
mkdir -p "$ND_DIR"
FOFN="$ND_DIR/input.fofn"
echo "$READS" > "$FOFN"

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

########################################
# 2. CANU
########################################
echo "Running Canu..."
Canu_DIR="$ROOT/canu"
mkdir -p "$Canu_DIR"
canu -p "$PROJECT_NAME" \
     -d "$Canu_DIR" \
     genomeSize="$GENOME_SIZE" \
     -nanopore-raw "$READS" \
     useGrid=false \
     minReadLength=500 \
     maxThreads=$THREADS 

########################################
# 3. FLYE
########################################
echo "Running Flye..."
export PATH="/home/silviad/minimap2-2.29_x64-linux:$PATH"

Flye_DIR="$ROOT/flye"
mkdir -p "$Flye_DIR"
flye --nano-corr "$READS" \
     --out-dir "$Flye_DIR" \
     --threads "$THREADS" \
     --genome-size "$GENOME_SIZE" \
     --min-overlap 1000 \
     --asm-coverage 50 | tee "$Flye_DIR/flye.log"
 

########################################
# 4. WTDBG2
########################################
echo "Running Wtdbg2..."
Wtdbg2_DIR="$ROOT/wtdbg2"
mkdir -p "$Wtdbg2_DIR"
TRIMMED_CANU_READS="${Canu_DIR}/${PROJECT_NAME}.trimmedReads.fasta.gz"
gunzip -c "$TRIMMED_CANU_READS" > "$Wtdbg2_DIR/corrected_reads.fasta"

wtdbg2 -i "$Wtdbg2_DIR/corrected_reads.fasta" -o "$Wtdbg2_DIR/wtdbg_assembly" -t $THREADS -g $GENOME_SIZE -p 0 -k 15 -AS 2 -s 0.05 -L 5000
wtpoa-cns -t $THREADS -i "$Wtdbg2_DIR/wtdbg_assembly.ctg.lay.gz" -fo "$Wtdbg2_DIR/wtdbg_assembly.fasta"

########################################
# 5. HIFIASM
########################################
echo "Running Hifiasm..."
Hifiasm_DIR="$ROOT/hifiasm"
mkdir -p "$Hifiasm_DIR"
hifiasm -o "$Hifiasm_DIR/$PROJECT_NAME" -t $THREADS --ont "$READS"

HIFIASM_GFA="$Hifiasm_DIR/${PROJECT_NAME}.bp.p_ctg.gfa"
HIFIASM_FASTA="$Hifiasm_DIR/${PROJECT_NAME}_from_gfa.fasta"
awk '/^S/{print ">"$2"\n"$3}' "$HIFIASM_GFA" > "$HIFIASM_FASTA"

echo "Deactivating assembly environment"
conda deactivate
echo "Activating polishing environment: $POLISHING_ENV"
conda activate "$POLISHING_ENV"

########################################
# 6. QUAST pre-polishing
########################################
echo "Running QUAST pre-polishing..."
QUAST_PRE="$ROOT/quast_results/pre_polishing"
mkdir -p "$QUAST_PRE"
declare -A assemblies=(
  ["canu"]="$Canu_DIR/${PROJECT_NAME}.contigs.fasta"
  ["nextdenovo"]="$(find "$ND_DIR" -name '*nd.asm*.fasta' | head -n1)"
  ["flye"]="$Flye_DIR/assembly.fasta"
  ["wtdbg2"]="$Wtdbg2_DIR/wtdbg_assembly.fasta"
  ["hifiasm"]="$HIFIASM_FASTA"
)

quast_cmd="quast.py -t $THREADS -o $QUAST_PRE"
for tool in "${!assemblies[@]}"; do
  if [ -s "${assemblies[$tool]}" ]; then
    quast_cmd+=" ${assemblies[$tool]}"
  else
    echo "Missing file for $tool: ${assemblies[$tool]}"
  fi
done
eval "$quast_cmd"

########################################
# 7. POLISHING: Racon + Medaka
########################################
echo "Running polishing with Racon and Medaka..."
POLISH_DIR="$ROOT/polished_assemblies"
mkdir -p "$POLISH_DIR"

for ASSEMBLY in "${assemblies[@]}"; do
  [ ! -s "$ASSEMBLY" ] && continue
  SAMPLE=$(basename "$ASSEMBLY" .fasta)
  OUT_DIR="$POLISH_DIR/$SAMPLE"
  mkdir -p "$OUT_DIR"
  
  echo "--> Polishing $SAMPLE"
  minimap2 -ax map-ont -t $THREADS "$ASSEMBLY" "$READS" > "$OUT_DIR/reads.sam"
  racon -t $THREADS "$READS" "$OUT_DIR/reads.sam" "$ASSEMBLY" > "$OUT_DIR/racon_polished.fasta"

  echo "Running Medaka env: $MEDAKA_ENV_PATH..."
  source "$MEDAKA_ENV_PATH/bin/activate"
  medaka_consensus -i "$READS" -d "$OUT_DIR/racon_polished.fasta" -o "$OUT_DIR/medaka" -m "$MODEL"
  deactivate

  [ ! -s "$OUT_DIR/medaka/consensus.fasta" ] && echo "Medaka did not generate output for $SAMPLE"
done

########################################
# 8. QUAST post-polishing
########################################
echo "Running QUAST post-polishing..."
QUAST_POST="$ROOT/quast_results/post_polishing"
mkdir -p "$QUAST_POST"
post_files=$(find "$POLISH_DIR" -name "consensus.fasta")

quast.py -t $THREADS -o "$QUAST_POST" $post_files

echo "Deactivating polishing environment"
conda deactivate


echo "Pipeline completed"
