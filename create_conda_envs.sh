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
conda create -y -n preprocessing_env -c bioconda -c conda-forge minimap2=2.28

echo "Activating preprocessing_env and installing pod5 via pip"
source activate preprocessing_env
pip install pod5

# Dorado install
echo "Downloading and extracting Dorado v0.9.6..."
echo "The current version of Dorado is suitable for both r9.4.1 and r10.4.1. Based on your sequencing run and configuration, you may need to use a newer version (v1.0.1)."
echo "Note: Make sure this version matches your ONT chemistry, if not, choose the correct version of Dorado on GitHub (https://github.com/nanoporetech/dorado/)."
read -p "Path to the directory to install Dorado: " DORADO_DIR
mkdir -p "$DORADO_DIR"
pushd "$DORADO_DIR"
curl -L "https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.6-linux-x64.tar.gz" -o dorado.tar.gz
tar -xzf dorado.tar.gz
rm dorado.tar.gz
./dorado-0.9.6-linux-x64/bin/dorado --version
popd
conda deactivate

DORADO_BIN="$DORADO_DIR/dorado-0.9.6-linux-x64/bin"
if [[ ":$PATH:" != *":$DORADO_BIN:"* ]]; then
  echo "Adding Dorado to PATH in ~/.bashrc"
  {
    echo ""
    echo "# Added by Dorado installer"
    echo "export PATH=\"\$PATH:$DORADO_BIN\""
  } >> ~/.bashrc
fi
source ~/.bashrc
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
