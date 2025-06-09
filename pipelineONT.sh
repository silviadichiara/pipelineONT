#!/bin/bash

read -p "Path to save the LOGfile (e.g. data/name/pipeline_ont.log): " LOGFILE
echo "=== RBH & Masking & Alignment Pipeline Started ===" | tee "$LOGFILE"
echo ""

# --- Preliminary suggestions ---
echo "Preliminary suggestions:" | tee -a "$LOGFILE"
echo " 1) Download reference gene FASTA files (nucleotide and protein) and GFF3 annotations from NCBI." | tee -a "$LOGFILE"
echo "    For example:" | tee -a "$LOGFILE"
echo "      - Visit https://www.ncbi.nlm.nih.gov/genome/ to download genomes of reference strains." | tee -a "$LOGFILE"
echo "      - Use wget for FTP downloads, e.g.:" | tee -a "$LOGFILE"
echo "        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/XXXX/YYYYY/GCF_XXXXX_genomic.fna.gz" | tee -a "$LOGFILE"
echo "      - Use Entrez Direct for batch downloads:" | tee -a "$LOGFILE"
echo "        https://www.ncbi.nlm.nih.gov/books/NBK179288/" | tee -a "$LOGFILE"
echo ""
echo " 2) Make sure all your scripts are in the same directory or update the paths below accordingly." | tee -a "$LOGFILE"
echo ""

read -p "Press ENTER to start the pipeline..."

run_step() {
    STEP_NAME=$1
    COMMAND=$2
    echo "$STEP_NAME" | tee -a "$LOGFILE"
    read -p "Do you want to run this step? (y/n): " answer
    if [[ "$answer" == "y" || "$answer" == "Y" ]]; then
        eval "$COMMAND 2>&1 | tee -a \"$LOGFILE\""
        echo ""
    else
        echo "Skipped: $STEP_NAME" | tee -a "$LOGFILE"
        echo ""
    fi
}

run_step "Step 1:  Dorado basecalling and trimming adapters..." "bash 1_dorado_bc_trim_int.sh"
run_step "Step 2: Taxonomic analysis and reads selection..." "bash 2_tax_sel_qc.sh"
run_step "Step 3: Assembly and polishing..." "bash 3_asm_pol.sh"
run_step "Step 4: Using Illumina reads for hybrid assembly and polishing..." "bash 4_hybrid.sh"
run_step "Step 5: Comparing denovo assembly to reference and performing variant calling..." "bash 5_comparation.sh"
run_step "Step 6: Performing Reciprocal Best Hit..." "bash 6_rbh.sh"
run_step "Step 7: Parse RBH glsearch output and create lists per identity..." "python3 7_rbh_lists.py"
run_step "Step 8: Extracting sequences based on TSV coordinates..." "bash 8_extraction_fastagenes.sh"
run_step "Step 9: Masking sequences with identity thresholds..." "bash 9_masking.sh"
run_step "Step 10: Extracting unaligned regions and aligning to whole NCBI dataset..." "bash extract_and_align_unmasked_regions.sh"
run_step "Step 11: Translating nucleotidic sequences to aminoacidic sequences..." "bash 11_nuctoaa.sh"

# === Final suggestions and pipeline end ===
echo "" | tee -a "$LOGFILE"
echo "=== Pipeline completed successfully ===" | tee -a "$LOGFILE"
echo "Next steps and tips:" | tee -a "$LOGFILE"
echo " - Use the Excel and TSV files generated for detailed analysis of identities and regions." | tee -a "$LOGFILE"
echo " - Download additional genomes or annotation files from:" | tee -a "$LOGFILE"
echo "     https://www.ncbi.nlm.nih.gov/genome/" | tee -a "$LOGFILE"
echo "     https://www.ncbi.nlm.nih.gov/assembly/" | tee -a "$LOGFILE"
echo " - Visualize your genomic data with tools like IGV" | tee -a "$LOGFILE"
echo ""
echo "Log file saved as: $LOGFILE"
echo "Happy analyzing!"
