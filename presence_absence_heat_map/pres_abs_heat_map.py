# --- Combined Script: Gene Presence/Absence CSV Creation and Heatmap Generation ---

import os
import pandas as pd
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# Usage:
# python combined_script.py --fasta_dir path/to/fasta_directory --gene_order_file path/to/gene_order.txt --species_file path/to/species_order.txt --output_csv path/to/output_file.csv --output_pdf path/to/output.pdf
#
# Arguments:
#   --fasta_dir: Directory containing gene FASTA files (.fasta or .faa)
#   --gene_order_file: Text file with desired gene order (one gene per line)
#   --species_file: Text file with desired species order (one species per line)
#   --output_csv: Output CSV file for gene presence/absence data
#   --output_pdf: Output PDF file for heatmap

def load_desired_gene_order(gene_order_file):
    """Load desired gene order from a file."""
    with open(gene_order_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def load_species_order(species_file):
    """Load species order from a file."""
    with open(species_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def create_presence_absence_csv(fasta_dir, gene_order_file, output_csv):
    """Process gene FASTA files and create a gene presence/absence CSV."""
    desired_gene_order = load_desired_gene_order(gene_order_file)
    data = {}
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fasta') or filename.endswith('.faa'):
            gene_name = os.path.splitext(filename)[0]
            print(f"Processing gene: {gene_name}")
            fasta_path = os.path.join(fasta_dir, filename)
            for record in SeqIO.parse(fasta_path, "fasta"):
                try:
                    parts = record.id.split('_')
                    if len(parts) > 4 and parts[4] == "rubra":
                        species_id = "_".join(parts[2:])
                    elif len(parts) > 1:
                        species_id = "_".join(parts[1:])
                    else:
                        print(f"Unexpected format in record ID: {record.id}")
                        continue
                    print(f"Identified species: {species_id} for gene: {gene_name}")
                    if species_id not in data:
                        data[species_id] = {}
                    data[species_id][gene_name] = 1
                except IndexError:
                    print(f"Error parsing species ID from record ID: {record.id}")
                    continue
    df = pd.DataFrame.from_dict(data, orient='index').fillna(0).astype(int)
    for gene in desired_gene_order:
        if gene not in df.columns:
            df[gene] = 0
    df = df[desired_gene_order]
    df.to_csv(output_csv)
    print(f"CSV file '{output_csv}' has been created.")
    return df

def plot_heatmap(input_csv, species_file, output_pdf):
    """Generate a heatmap of gene presence/absence data."""
    df = pd.read_csv(input_csv, index_col=0)
    species_order = load_species_order(species_file)
    for species in species_order:
        if species not in df.index:
            df.loc[species] = 0
    df = df.reindex(species_order)
    plt.figure(figsize=(10, 15))
    heatmap = sns.heatmap(df, cmap="YlGnBu", linewidths=0.5, linecolor='gray', cbar_kws={"shrink": .5})
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, fontsize=8)
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=90, fontsize=8)
    plt.tight_layout()
    plt.savefig(output_pdf, format='pdf')
    print(f"Heatmap saved as '{output_pdf}'.")
    plt.show()
    return df

def main():
    parser = argparse.ArgumentParser(description='Create gene presence/absence CSV and generate heatmap.')
    parser.add_argument('--fasta_dir', type=str, required=True, help='Directory containing gene FASTA files.')
    parser.add_argument('--gene_order_file', type=str, required=True, help='File with desired gene order.')
    parser.add_argument('--species_file', type=str, required=True, help='File with desired species order.')
    parser.add_argument('--output_csv', type=str, required=True, help='Output CSV file.')
    parser.add_argument('--output_pdf', type=str, required=True, help='Output PDF file for heatmap.')
    args = parser.parse_args()

    # Step 1: Create presence/absence CSV
    df_csv = create_presence_absence_csv(args.fasta_dir, args.gene_order_file, args.output_csv)
    print("\nPresence/Absence DataFrame:")
    print(df_csv)

    # Step 2: Generate heatmap
    df_heatmap = plot_heatmap(args.output_csv, args.species_file, args.output_pdf)
    print("\nHeatmap DataFrame (after species reordering):")
    print(df_heatmap)

if __name__ == "__main__":
    main()
