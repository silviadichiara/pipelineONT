#!/usr/bin/env bash
conda activate fasta3_env
# Controllo che fasty36 sia installato
if ! command -v fasty36 &> /dev/null; then
  echo "Error: fasty36 command not found. Please install it before running this script."
  exit 1
fi

# Input da utente
read -p "Enter path to query nucleotide file (e.g., denovo_genes_seq.fasta): " QUERY
read -p "Enter path to reference protein file (e.g., vr1814_geni_aa.fasta): " REF_PROTEIN
read -p "Enter path to readable output file (e.g., fasty_output_readable.txt): " OUTFMT0
read -p "Enter path to tabular output file (e.g., fasty_output_tabular.tsv): " OUTFMT6
read -p "Enter path to summary output file (e.g., fasty_best_hits_summary.tsv): " SUMMARY

echo -e "\n[1] Running fasty36 in readable format..."
fasty36 -E 1e-5 -m 0 "$QUERY" "$REF_PROTEIN" > "$OUTFMT0" || {
    echo "Error running fasty36 (readable output)"; exit 1;
}

echo -e "\n[2] Running fasty36 in tabular format..."
fasty36 -E 1e-5 -m 8 "$QUERY" "$REF_PROTEIN" > "$OUTFMT6" || {
    echo "Error running fasty36 (tabular output)"; exit 1;
}

echo -e "\n[3] Extracting best hits from tabular output..."
awk '
BEGIN {
    FS = OFS = "\t"
    print "query_id", "subject_id", "percent_identity", "alignment_length", "mismatches", "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit_score"
}
!/^#/ {
    if (!seen[$1]++) print $0
}
' "$OUTFMT6" > "$SUMMARY"

echo -e "\nSummary completed:"
echo " - Readable format saved to: $OUTFMT0"
echo " - Tabular format saved to: $OUTFMT6"
echo " - Best hits summary saved to: $SUMMARY"
conda deactivate
