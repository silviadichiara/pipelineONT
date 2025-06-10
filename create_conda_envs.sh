#!/bin/bash

set -e  
set -o pipefail

if ! command -v conda &> /dev/null; then
  echo "Error: Conda is not installed or not in PATH."
  exit 1
fi

eval "$(conda shell.bash hook)"

echo "Starting Conda environment setup..."

# === 1. preprocessing_env ===
echo "Creating environment: preprocessing_env"
conda create -y -n preprocessing_env -c bioconda -c conda-forge \
  minimap2=2.28

echo "Activating preprocessing_env and installing pod5 via pip"
conda activate preprocessing_env
pip install pod5

echo "Downloading and extracting Dorado..."
read -p "Path to the directory to install Dorado: " DORADO_DIR
mkdir -p $DORADO_DIR
cd $DORADO_DIR
curl -L "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.0.1-linux-x64.tar.gz" -o dorado.tar.gz
tar -xzf dorado.tar.gz
rm dorado.tar.gz
./dorado-1.0.1-linux-x64/bin/dorado --version
cd -
conda deactivate

# === 2. kraken2_env ===
echo "Creating environment: kraken2_env"
conda create -y -n kraken2_env -c bioconda -c conda-forge kraken2 seqtk

# === 3. assembly_env ===
echo "Creating environment: assembly_env"
conda create -y -n assembly_env -c bioconda -c conda-forge \
  flye canu hifiasm wtdbg paralleltask nextdenovo

# === 4. polishing_env ===
echo "Creating environment: polishing_env"
conda create -y -n polishing_env -c bioconda -c nanoporetech -c conda-forge \
  medaka racon quast

# === 5. qc_env ===
echo "Creating environment: qc_env"
conda create -y -n qc_env -c bioconda -c conda-forge \
  nanoplot fastqc qualimap bcftools samtools

# === 6. illuminareads_env ===
echo "Creating environment: illuminareads_env"
conda create -y -n illuminareads_env -c bioconda -c conda-forge \
  bwa pilon freebayes

# === 7. mummer_env ===
echo "Creating environment: mummer_env"
conda create -y -n mummer_env -c bioconda -c conda-forge \
  mummer bedtools blast

# === 8. fasta3_env ===
echo "Creating environment: fasta3_env"
conda create -y -n fasta3_env -c conda-forge -c bioconda \
  pandas openpyxl biopython fasta3

echo "All environments have been successfully created!"
