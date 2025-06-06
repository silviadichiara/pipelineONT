#!/bin/bash

read -p "Enter path to the reference FASTA file: " FASTA
read -p "Enter path to the TSV file (RBH_all.tsv): " TSV
read -p "Enter desired output FASTA file path: " OUTPUT

# Clear output file
> "$OUTPUT"

awk -v fasta="$FASTA" -v tsv="$TSV" -v out="$OUTPUT" '
BEGIN {
    FS="\t"
    seq_len = 0
    seq = ""
    header_found = 0

    # Load the first sequence from FASTA file
    while ((getline line < fasta) > 0) {
        if (line ~ /^>/) {
            if (header_found) break;  # stop if second header found (only first seq loaded)
            header_found = 1
            continue
        }
        seq = seq line
    }
    close(fasta)
    seq_len = length(seq)
    if (seq_len == 0) {
        print "Error: No sequence loaded from FASTA."
        exit 1
    }
}

FNR > 1 {
    id = $1
    start = $(NF-1)
    end = $NF

    if (start !~ /^[0-9]+$/ || end !~ /^[0-9]+$/) {
        print "Warning: invalid coordinates for " id ", skipping"
        next
    }

    if (start > end) {
        tmp = start
        start = end
        end = tmp
    }

    if (start < 1 || end > seq_len) {
        print "Warning: coordinates out of bounds for " id ", skipping"
        next
    }

    len = end - start + 1
    fragment = substr(seq, start, len)

    printf(">%s %d-%d len=%d\n%s\n", id, start, end, len, fragment) >> out
    close(out)  # flush output after each write
}
' "$TSV"

echo "Sequences successfully extracted to: $OUTPUT"
