#!/bin/bash

set -e
set -o pipefail

# Check if conda is available
if ! command -v conda &> /dev/null; then
  echo "Error: Conda is not installed or not in PATH."
  exit 1
fi

# Initialize conda in current shell
eval "$(conda shell.bash hook)"
echo "Starting Conda environment setup..."

# === 1. preprocessing_env ===
echo "Creating environment: preprocessing_env"
conda create -y -n preprocessing_env -c bioconda -c conda-forge minimap2=2.28

echo "Installing pod5 in preprocessing_env"
conda activate preprocessing_env
pip install pod5
conda deactivate

# === Dorado install ===
echo "Downloading and extracting Dorado v0.9.6..."
echo "Note: For different ONT chemistries, see https://github.com/nanoporetech/dorado"
read -p "Path to the directory where Dorado should be installed: " DORADO_DIR
mkdir -p "$DORADO_DIR"
pushd "$DORADO_DIR"
curl -L "https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.6-linux-x64.tar.gz" -o dorado.tar.gz
tar -xzf dorado.tar.gz
rm dorado.tar.gz
./dorado-0.9.6-linux-x64/bin/dorado --version
popd

# Add Dorado to PATH if not already added
DORADO_BIN="$DORADO_DIR/dorado-0.9.6-linux-x64/bin"
if [[ ":$PATH:" != *":$DORADO_BIN:"* ]]; then
  echo "Adding Dorado to PATH in ~/.bashrc"
  {
    echo ""
    echo "# Added by Dorado installer"
    echo "export PATH=\"\$PATH:$DORADO_BIN\""
  } >> ~/.bashrc
fi
export PATH="$PATH:$DORADO_BIN"

# === 2. kraken2_env ===
echo "Creating environment: kraken2_env"
conda create -y -n kraken2_env -c bioconda -c conda-forge kraken2 seqtk

# === 3. assembly_env ===
echo "Creating environment: assembly_env"
conda create -y -n assembly_env -c bioconda -c conda-forge \
  flye hifiasm wtdbg paralleltask

conda activate assembly_env
conda install -c conda-forge openjdk=8

echo "Downloading and extracting Canu v2.3..."
read -p "Path to the directory where Canu should be installed: " CANU_DIR

mkdir -p "$CANU_DIR"
pushd "$CANU_DIR"
curl -LRO https://github.com/marbl/canu/releases/download/v2.3/canu-2.3.Linux-amd64.tar.xz
tar -xJf canu-2.3.Linux-amd64.tar.xz
rm canu-2.3.Linux-amd64.tar.xz

./canu-2.3/Linux-amd64/bin/canu --version

popd

CANU_BIN="$CANU_DIR/canu-2.3/Linux-amd64/bin"
if [[ ":$PATH:" != *":$CANU_BIN:"* ]]; then
  echo "Adding Canu to PATH in ~/.bashrc"
  {
    echo ""
    echo "# Added by Canu installer"
    echo "export PATH=\"\$PATH:$CANU_BIN\""
  } >> ~/.bashrc
fi

export PATH="$PATH:$CANU_BIN"


# === 4. nextdenovo in separate env ===
echo "Creating environment: ndn_env (NextDenovo)"
conda create -y -n ndn_env -c bioconda nextdenovo

# === 5. polishing_env ===
echo "Creating environment: polishing_env"
conda create -y -n polishing_env -c bioconda -c nanoporetech -c conda-forge \
  medaka racon quast

# === 6. qc_env ===
echo "Creating environment: qc_env"
conda create -y -n qc_env -c bioconda -c conda-forge \
  nanoplot fastqc qualimap bcftools samtools

# === 7. illuminareads_env ===
echo "Creating environment: illuminareads_env"
conda create -y -n illuminareads_env -c bioconda -c conda-forge \
  bwa pilon freebayes

# === 8. mummer_env ===
echo "Creating environment: mummer_env"
conda create -y -n mummer_env -c bioconda -c conda-forge \
  mummer bedtools blast

# === 9. fasta3_env ===
echo "Creating environment: fasta3_env"
conda create -y -n fasta3_env -c conda-forge -c bioconda \
  pandas openpyxl biopython fasta3

echo "All environments have been successfully created."
