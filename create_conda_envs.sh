#!/bin/bash

set -e  
set -o pipefail

if ! command -v conda &> /dev/null; then
  echo "Error: Conda is not installed or not in PATH."
  exit 1
fi

eval "$(conda shell.bash hook)"

echo "Starting Conda environment setup..."

# 1. preprocessing_env 
echo "Creating environment: preprocessing_env"
conda create -y -n preprocessing_env -c bioconda -c conda-forge \
  minimap2=2.28 dorado python=3.10
conda activate preprocessing_env
echo "Installing pod5 via pip..."
pip install pod5
conda deactivate 

# 2. kraken2_env
conda create -y -n kraken2_env -c bioconda -c conda-forge kraken2 seqtk

# 3. assembly_env
conda create -y -n assembly_env -c bioconda -c conda-forge \
  flye canu hifiasm wtdbg paralleltask nextdenovo

# 4. polishing_env
conda create -y -n polishing_env -c bioconda -c nanoporetech -c conda-forge \
  medaka racon quast

# 5. qc_env
conda create -y -n qc_env -c bioconda -c conda-forge \
  nanoplot fastqc qualimap bcftools samtools

# 6. illuminareads_env
conda create -y -n illuminareads_env -c bioconda -c conda-forge \
  bwa pilon freebayes

# 7. mummer_env
conda create -y -n mummer_env -c bioconda -c conda-forge \
  mummer bedtools blast
  
# 8. fasta3_env
conda install -y -c conda-forge -c bioconda \
  pandas openpyxl biopython fasta3

echo "All environments have been successfully created!"
