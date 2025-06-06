# -*- coding: utf-8 -*-

import os
import re
import pandas as pd

# ==== Input paths from user ====
contig_vs_genes_dir = input("Enter path to 'contig_vs_genes' directory: ").strip()
genes_vs_contig_dir = input("Enter path to 'genes_vs_contig' directory: ").strip()
output_dir = input("Enter path to output directory: ").strip()
os.makedirs(output_dir, exist_ok=True)

# ==== Function to parse glsearch output ====
def parse_glsearch(filepath):
    identity_pattern = re.compile(r'([\d\.]+)% identity')
    coords_pattern = re.compile(r'overlap \(\d+-\d+:(\d+)-(\d+)\)')
    hit_name_pattern = re.compile(r'^>>(\S+)')

    identity = None
    coords = None
    hit_name = None

    try:
        with open(filepath, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if hit_name is None and line.startswith(">>"):
                    hit_match = hit_name_pattern.match(line)
                    if hit_match:
                        hit_name = hit_match.group(1)
                if identity is None:
                    id_match = identity_pattern.search(line)
                    if id_match:
                        identity = float(id_match.group(1))
                if coords is None:
                    coord_match = coords_pattern.search(line)
                    if coord_match:
                        start = int(coord_match.group(1))
                        end = int(coord_match.group(2))
                        coords = (min(start, end), max(start, end))
                if identity is not None and coords is not None and hit_name is not None:
                    break
    except Exception as e:
        print(f"Error reading {filepath}: {e}")

    if identity is None or coords is None or hit_name is None:
        print(f"Warning: Missing data in {filepath} (identity={identity}, coords={coords}, hit={hit_name})")

    return identity, coords, hit_name

# ==== Helper functions for RBH determination ====
def overlaps(coords1, coords2, tol=5):
    """Check if two coordinate intervals overlap within a tolerance."""
    return not (coords1[1] + tol < coords2[0] or coords2[1] + tol < coords1[0])

def is_rbh(c_identity, g_identity, c_coords, g_coords, identity_tol=1.0, coord_tol=5):
    """Determine if two hits form a Reciprocal Best Hit (RBH) with tolerances."""
    if None in (c_identity, g_identity, c_coords, g_coords):
        return False
    identity_close = abs(c_identity - g_identity) <= identity_tol
    coords_overlap = overlaps(c_coords, g_coords, tol=coord_tol)
    return identity_close and coords_overlap

# ==== Process all matching result files ====
rbh_records = []
total_files = 0
skipped_files = 0

for filename in os.listdir(contig_vs_genes_dir):
    if not filename.endswith(".glsearch.out"):
        continue

    total_files += 1

    contig_file = os.path.join(contig_vs_genes_dir, filename)
    gene_file = os.path.join(genes_vs_contig_dir, filename)

    if not os.path.exists(gene_file):
        print(f"Warning: Matching gene file not found for {filename}, skipping.")
        skipped_files += 1
        continue

    # Extract base name removing extension
    basename = filename.replace(".glsearch.out", "")

    contig_identity, contig_coords, contig_hit = parse_glsearch(contig_file)
    gene_identity, gene_coords, gene_hit = parse_glsearch(gene_file)

    if is_rbh(contig_identity, gene_identity, contig_coords, gene_coords):
        avg_identity = (contig_identity + gene_identity) / 2
        rbh_records.append({
            "Gene": basename,
            "Identity": avg_identity,
            "Contig_hit": contig_hit,
            "Contig_start": contig_coords[0],
            "Contig_end": contig_coords[1],
            "Gene_hit": gene_hit,
            "Gene_start": gene_coords[0],
            "Gene_end": gene_coords[1]
        })
    else:
        print(f"Skipped non-RBH: {basename} "
              f"(Identities: {contig_identity}, {gene_identity}; "
              f"Coords: {contig_coords}, {gene_coords})")

# ==== Create DataFrame and filter ====
df = pd.DataFrame(rbh_records)

df_100 = df[df["Identity"] >= 99.99]  # Approximating 100% due to float precision
df_95_99 = df[(df["Identity"] >= 95.0) & (df["Identity"] < 99.99)]
df_below_95 = df[df["Identity"] < 95.0]

# ==== Write Excel file ====
excel_path = os.path.join(output_dir, "RBH_results.xlsx")
with pd.ExcelWriter(excel_path) as writer:
    df.to_excel(writer, sheet_name="All_RBH", index=False)
    df_100.to_excel(writer, sheet_name="Identity_100", index=False)
    df_95_99.to_excel(writer, sheet_name="Identity_95_99", index=False)
    df_below_95.to_excel(writer, sheet_name="Identity_below_95", index=False)

# ==== Write TSV files ====
df.to_csv(os.path.join(output_dir, "RBH_all.tsv"), sep="\t", index=False)
df_100.to_csv(os.path.join(output_dir, "RBH_100.tsv"), sep="\t", index=False)
df_95_99.to_csv(os.path.join(output_dir, "RBH_95_99.tsv"), sep="\t", index=False)
df_below_95.to_csv(os.path.join(output_dir, "RBH_below_95.tsv"), sep="\t", index=False)

# ==== Print summary ====
print(f"Processed {total_files} files.")
print(f"Skipped {skipped_files} files due to missing pairs.")
print(f"Total RBHs found: {len(df)}")
print(f" - 100% identity (≥99.99): {len(df_100)}")
print(f" - 95-99% identity: {len(df_95_99)}")
print(f" - Below 95% identity: {len(df_below_95)}")
print(f"Excel file saved to: {excel_path}")
print("TSV files saved in directory:", output_dir)
