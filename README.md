## 🧬Hybrid Comparative Genomics Pipeline

This repository provides an automated pipeline for comparative genomic analysis using a hybrid approach integrating:
- Oxford Nanopore Technologies (ONT) long-read DNA sequencing
- Illumina RNA-seq short-read transcriptome data

It streamlines analysis from raw signal to polished genome assemblies, including comparison and evaluation metrics.

## 💻 Features
This pipeline includes scripts and configurations for:

- **Conversion**: fast5 → pod5 using pod5tools

- **Basecalling**: performed using Dorado

- **Taxonomic analysis**: Kraken2

- **Quality Control**: Nanoplot, FastQC, Qualimap

- **Trimming & Filtering**: fastp, seqtk

- **De novo Assembly**: Flye, Canu, wtdbg2, NextDenovo, Hifiasm
  
- **Comparison & Evaluation**: QUAST, MUMmer, Minimap2, BLAST, GLSEARCH, BCFtools
    
- **Polishing**:

  -*DNA-based*: Racon, Medaka

  -*RNA-based*: BWA, Pilon


## ✅ Requirements
- Bash (Linux environment)
- Conda or Miniconda installed
- Compatible bioinformatics tools installed and accessible in your $PATH (see environment setup below)

## 🛠 Installation and Setup
**Create the Conda environments and install tools**:
Run the script create_conda_envs.sh

## 🚀 Running the Pipeline
Run the main orchestration script
The pipeline is controlled by a master bash script which runs all analysis steps and logs progress:
bash pipeline_ONT.sh

## 📊 Outputs
- Polished genome assemblies and their stats and reports;
- Alignment to reference genes and RBH summary Excel and TSV files;
- Nucleotidic and aminoacidic sequences from denovo assembly FASTA file;
- Alignment logs and summary reports.

## 📚 Additional Notes
For reference genomes and annotation data, download FASTA and GFF3 files from NCBI:
- NCBI Genomes
- Use tools like wget or Entrez Direct for batch downloads.
- Visualize your results with genome browsers like IGV or Artemis.

## 🧠 Citation
If you use this pipeline in your research, please cite the relevant tools listed above and credit this repository.

## 📫 Contact
For questions or contributions, please open an issue or submit a pull request.
