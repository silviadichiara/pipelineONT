#!/bin/bash
conda activate fasta3_env
read -p "Enter path to TSV file: " TSV
read -p "Enter path to FASTA file: " FASTA

PREFIX=$(basename "$FASTA")
PREFIX=${PREFIX%%.*}
OUTDIR=$(dirname "$FASTA")/masked_output
mkdir -p "$OUTDIR"

python3 - <<END
import csv
from decimal import Decimal, InvalidOperation
import os

tsv_file = "${TSV}"
fasta_file = "${FASTA}"
prefix = "${PREFIX}"
output_dir = "${OUTDIR}"

def parse_tsv(tsv_path, min_id=Decimal("0.0"), max_id=Decimal("100.0")):
    regions = []
    with open(tsv_path, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            try:
                identity = Decimal(row['Identity'].strip())
                start = int(row['Contig_start'])
                end = int(row['Contig_end'])
            except (InvalidOperation, ValueError, KeyError):
                continue
            if min_id <= identity <= max_id:
                regions.append((start, end))
    return regions

def mask_sequence(sequence, regions):
    seq_list = list(sequence)
    for start, end in regions:
        for i in range(start - 1, end):
            seq_list[i] = 'N'
    return ''.join(seq_list)

def write_fasta(header, sequence, output_path):
    with open(output_path, 'w') as out:
        out.write(f"{header}\n")
        for i in range(0, len(sequence), 80):
            out.write(sequence[i:i+80] + '\n')

def load_fasta(fasta_path):
    with open(fasta_path, 'r') as f:
        header = f.readline().strip()
        seq_lines = []
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                break
            seq_lines.append(line)
    return header, ''.join(seq_lines)

header, sequence = load_fasta(fasta_file)

regions_100 = parse_tsv(tsv_file, Decimal("100.0"), Decimal("100.0"))
masked_100 = mask_sequence(sequence, regions_100)
write_fasta(header, masked_100, os.path.join(output_dir, f"{prefix}_masked_100.fasta"))

regions_95_999 = parse_tsv(tsv_file, Decimal("95.0"), Decimal("99.9"))
masked_95_999 = mask_sequence(sequence, regions_95_999)
write_fasta(header, masked_95_999, os.path.join(output_dir, f"{prefix}_masked_95_99.9.fasta"))

regions_all = parse_tsv(tsv_file, Decimal("0.0"), Decimal("100.0"))
masked_all = mask_sequence(sequence, regions_all)
write_fasta(header, masked_all, os.path.join(output_dir, f"{prefix}_masked_all.fasta"))
conda deactivate
print("Masking completed. Outputs saved in:", output_dir)
END
