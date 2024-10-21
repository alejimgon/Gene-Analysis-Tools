import os
import pandas as pd
from Bio import SeqIO
import argparse

# Use: python presence_absence.py path/to/fasta_directory path/to/gene_order.txt path/to/output_file.csv
# Input files:
# Input FASTA Files: Place your gene FASTA files in a directory. Each file should have a .fasta or .faa extension.
# Gene Order File: Create a text file named gene_order.txt. Each line should contain one gene that you want in the
# specified order.

def load_desired_gene_order(gene_order_file):
    """Load desired gene order from a file."""
    with open(gene_order_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def main(fasta_dir, gene_order_file, output_csv):
    # Load the desired gene order from the specified file
    desired_gene_order = load_desired_gene_order(gene_order_file)

    # Initialize an empty dictionary to store the presence/absence data
    data = {}

    # Loop over each FASTA file in the directory
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fasta') or filename.endswith('.faa'):
            gene_name = os.path.splitext(filename)[0]
            print(f"Processing gene: {gene_name}")

            # Parse the FASTA file
            fasta_path = os.path.join(fasta_dir, filename)
            for record in SeqIO.parse(fasta_path, "fasta"):
                try:
                    parts = record.id.split('_')
                    if len(parts) > 4 and parts[4] == "rubra":
                        species_id = "_".join(parts[2:])
                    elif len(parts) > 1:
                        species_id = "_".join(parts[1:])  # Join all parts after the first one
                    else:
                        print(f"Unexpected format in record ID: {record.id}")
                        continue
                    print(f"Identified species: {species_id} for gene: {gene_name}")

                    if species_id not in data:
                        data[species_id] = {}
                    data[species_id][gene_name] = 1  # Mark gene as present

                except IndexError:
                    print(f"Error parsing species ID from record ID: {record.id}")
                    continue

    # Convert the data dictionary to a DataFrame
    df = pd.DataFrame.from_dict(data, orient='index').fillna(0).astype(int)

    # Reorder the columns based on the desired gene order
    # Ensure all desired genes are included, even if they don't appear in the data
    for gene in desired_gene_order:
        if gene not in df.columns:
            df[gene] = 0  # Add missing genes as columns with 0s

    # Reorder the columns
    df = df[desired_gene_order]

    # Save the DataFrame to a CSV file
    df.to_csv(output_csv)

    print(f"CSV file '{output_csv}' has been created.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process gene FASTA files and create a gene presence/absence CSV.')
    parser.add_argument('fasta_dir', type=str, help='Directory containing the gene FASTA files.')
    parser.add_argument('gene_order_file', type=str, help='File containing the desired gene order (one gene per line).')
    parser.add_argument('output_csv', type=str, help='Output path for the CSV file.')

    args = parser.parse_args()

    main(args.fasta_dir, args.gene_order_file, args.output_csv)