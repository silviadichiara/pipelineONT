## 🧬 Hybrid Comparative Genomics Pipeline

This repository provides an automated pipeline for comparative genomic analysis using a hybrid approach that integrates:
Oxford Nanopore Technologies (ONT) long-read DNA sequencing
Illumina RNA-seq short-read transcriptome data
It streamlines the analysis from raw signal to polished genome assemblies, incorporating comparison and evaluation metrics.

## 💻 Features

This pipeline includes scripts and configurations for:

**Conversion**: fast5 → pod5 using pod5tools

**Basecalling**: performed using Dorado

**Taxonomic analysis**: Kraken2

**Quality Control**: Nanoplot, FastQC, Qualimap

**Trimming & Filtering**: fastp, seqtk

**De novo Assembly**: Flye, Canu, wtdbg2, NextDenovo, Hifiasm

**Polishing**:
- **DNA**-based: Racon, Medaka
- **RNA**-based: BWA, Pilon
  
**Comparison & Evaluation**: QUAST, MUMmer, Minimap2, BLAST, GLSEARCH, BCFtools

## ✅ Requirements
- Bash (Linux environment)
- Compatible tools installed and accessible in $PATH

## 📊 Outputs
- Raw and filtered reads
- Multiple polished assemblies
- QUAST reports (pre- and post-polishing)
- Comparative alignment summaries
- Final variant calls and annotations

## 🧠 Citation
If you use this pipeline in your research, please cite the tools listed above and credit this repository.

## 📫 Contact
Created and maintained by Silvia Di Chiara
For questions or contributions, please open an issue or submit a pull request.
