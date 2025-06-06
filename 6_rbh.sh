#!/bin/bash

read -p "Enter path to reference genes FASTA file: " REFERENCE_GENES
read -p "Enter path to de novo genome FASTA file: " DENOVO_GENOME
read -p "Enter path to output directory: " OUTPUT_DIR

mkdir -p "$OUTPUT_DIR/contig_vs_genes"
mkdir -p "$OUTPUT_DIR/genes_vs_contig"

if ! command -v glsearch36 &> /dev/null; then
    echo "Error: glsearch36 not found. Please install the FASTA suite."
    exit 1
fi

write_temp_fasta() {
    local header="$1"
    local seq="$2"
    local file="$3"
    echo -e ">$header\n$seq" > "$file"
}

sanitize_filename() {
    local name="$1"
    # Sostituisce caratteri problematici con underscore
    echo "$name" | tr ' /|[]()=' '______' | tr -cd '[:alnum:]_.-'
}

# 1. Allineamento: contig (query) contro ogni gene (database)
align_contig_to_genes() {
    local CONTIG_FASTA=$1
    local GENES_FASTA=$2
    local OUT_DIR=$3

    # Leggo il contig solo una volta (assumo un unico record)
    local CONTIG_SEQ=$(sed -n '2,$p' "$CONTIG_FASTA" | tr -d '\n')

    local HEADER=""
    local SEQ=""

    while read -r line; do
        if [[ $line == ">"* ]]; then
            if [[ -n "$HEADER" && -n "$SEQ" ]]; then
                local SAFE_HEADER=$(sanitize_filename "$HEADER")
                write_temp_fasta "$HEADER" "$SEQ" temp_gene.fasta
                # Qui il gene è database, contig è query
                glsearch36 -m 0 -n temp_gene.fasta <(echo -e ">contig\n$CONTIG_SEQ") > "$OUT_DIR/${SAFE_HEADER}.glsearch.out" || {
                    echo "Error during contig_vs_genes for gene $HEADER"
                    rm -f temp_gene.fasta
                    exit 1
                }
            fi
            HEADER="${line#>}"
            SEQ=""
        else
            SEQ+="$line"
        fi
    done < "$GENES_FASTA"

    # Ultimo record
    if [[ -n "$HEADER" && -n "$SEQ" ]]; then
        local SAFE_HEADER=$(sanitize_filename "$HEADER")
        write_temp_fasta "$HEADER" "$SEQ" temp_gene.fasta
        glsearch36 -m 0 -n temp_gene.fasta <(echo -e ">contig\n$CONTIG_SEQ") > "$OUT_DIR/${SAFE_HEADER}.glsearch.out" || {
            echo "Error during contig_vs_genes for gene $HEADER"
            rm -f temp_gene.fasta
            exit 1
        }
    fi

    rm -f temp_gene.fasta
}

# 2. Allineamento: ogni gene (query) contro contig (database)
align_genes_to_contig() {
    local GENES_FASTA=$1
    local CONTIG_FASTA=$2
    local OUT_DIR=$3

    # Leggo il contig solo una volta (assumo un unico record)
    local CONTIG_SEQ=$(sed -n '2,$p' "$CONTIG_FASTA" | tr -d '\n')
    write_temp_fasta "contig" "$CONTIG_SEQ" temp_contig.fasta

    local HEADER=""
    local SEQ=""

    while read -r line; do
        if [[ $line == ">"* ]]; then
            if [[ -n "$HEADER" && -n "$SEQ" ]]; then
                local SAFE_HEADER=$(sanitize_filename "$HEADER")
                write_temp_fasta "$HEADER" "$SEQ" temp_gene.fasta
                # Qui il contig è database, gene è query
                glsearch36 -m 0 -n temp_contig.fasta temp_gene.fasta > "$OUT_DIR/${SAFE_HEADER}.glsearch.out" || {
                    echo "Error during genes_vs_contig for gene $HEADER"
                    rm -f temp_gene.fasta temp_contig.fasta
                    exit 1
                }
            fi
            HEADER="${line#>}"
            SEQ=""
        else
            SEQ+="$line"
        fi
    done < "$GENES_FASTA"

    # Ultimo record
    if [[ -n "$HEADER" && -n "$SEQ" ]]; then
        local SAFE_HEADER=$(sanitize_filename "$HEADER")
        write_temp_fasta "$HEADER" "$SEQ" temp_gene.fasta
        glsearch36 -m 0 -n temp_contig.fasta temp_gene.fasta > "$OUT_DIR/${SAFE_HEADER}.glsearch.out" || {
            echo "Error during genes_vs_contig for gene $HEADER"
            rm -f temp_gene.fasta temp_contig.fasta
            exit 1
        }
    fi

    rm -f temp_gene.fasta temp_contig.fasta
}

# Esecuzione
echo "Running contig_vs_genes alignment (contig query vs each gene database)..."
align_contig_to_genes "$DENOVO_GENOME" "$REFERENCE_GENES" "$OUTPUT_DIR/contig_vs_genes"

echo "Running genes_vs_contig alignment (each gene query vs contig database)..."
align_genes_to_contig "$REFERENCE_GENES" "$DENOVO_GENOME" "$OUTPUT_DIR/genes_vs_contig"

echo "Alignments completed."
echo "Results saved in:"
echo " - $OUTPUT_DIR/contig_vs_genes"
echo " - $OUTPUT_DIR/genes_vs_contig"
