#!/bin/bash

set -e  
set -o pipefail

if ! command -v conda &> /dev/null; then
  error "Conda is not installed or not in PATH."
  exit 1
fi

log "Starting Conda environment setup..."

# === 1. preprocessing_env ===
log "Creating environment: preprocessing_env"
conda create -y -n preprocessing_env -c bioconda -c conda-forge -c hcc \
  lib-pod5 minimap2=2.28 dorado

# === 2. kraken2_env ===
log "Creating environment: kraken2_env"
conda create -y -n kraken2_env -c bioconda -c conda-forge \
  kraken2 seqtk

# === 3. assembly_env ===
log "Creating environment: assembly_env"
conda create -y -n assembly_env -c bioconda -c conda-forge \
  flye canu hifiasm wtdbg paralleltask nextdenovo

# === 4. polishing_env ===
log "Creating environment: polishing_env"
conda create -y -n polishing_env -c bioconda -c nanoporetech -c conda-forge \
  medaka racon quast

# === 5. qc_env ===
log "Creating environment: qc_env"
conda create -y -n qc_env -c bioconda -c conda-forge \
  nanoplot fastqc qualimap bcftools samtools

# === 6. illuminareads_env ===
log "Creating environment: illuminareads_env"
conda create -y -n illuminareads_env -c bioconda -c conda-forge \
  bwa pilon freebayes

# === 7. mummer_env ===
log "Creating environment: mummer_env"
conda create -y -n mummer_env -c bioconda -c conda-forge \
  mummer bedtools blast

conda install -c conda-forge -c bioconda \
  pandas openpyxl biopython fasta3

log "All environments have been successfully created!"
