import csv
import json
import os
import argparse
from Bio import SeqIO

# Use: python interproscan_pfam_filtering.py /path/to/your/folder path/to/csv_filename
# Inputs:
# JSON Files: Each protein should have a corresponding JSON file named <protein_name>.json.
# FASTA Files: Each protein should have a corresponding FASTA file named <protein_name>.faa.
# CSV File: Create a CSV file named protein_pfam.csv (or any name you choose) that contains the protein names in the
# first column and their respective PFAM domains in subsequent columns.

# Function to extract PFAM domains from each result in JSON
def extract_pfam_domains(result):
    """Extract PFAM domains from a result in JSON data."""
    pfam_domains = set()
    for match in result.get('matches', []):
        signature = match.get('signature', {})
        library = signature.get('signatureLibraryRelease', {}).get('library')
        if library == 'PFAM':
            pfam_domains.add(signature.get('accession'))
    return pfam_domains


def main(directory, csv_filename):
    # Read the CSV file and prepare a dictionary with protein names and their PFAM domains
    protein_pfam_dict = {}
    with open(os.path.join(directory, csv_filename), 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            protein_name = row[0]
            pfam_domains = set(row[1:])
            if pfam_domains:
                protein_pfam_dict[protein_name] = pfam_domains
            else:
                print(f"No PFAM domains found for {protein_name} in the CSV file.")

    # Prepare a dictionary to store sequence IDs for each protein based on PFAM domains
    protein_sequence_ids = {protein: set() for protein in protein_pfam_dict}

    # Process each protein
    for protein, pfam_domains_in_csv in protein_pfam_dict.items():
        # Construct the file paths
        json_file_path = os.path.join(directory, f'{protein}.json')
        fasta_file_path = os.path.join(directory, f'{protein}.faa')
        output_fasta_path = os.path.join(directory, f'{protein}_filtered.faa')

        if not os.path.exists(json_file_path):
            print(f"JSON file not found for {protein}, skipping.")
            continue

        if not os.path.exists(fasta_file_path):
            print(f"FASTA file not found for {protein}, skipping.")
            continue

        # Read the JSON file
        try:
            with open(json_file_path, 'r') as file:
                data = json.load(file)
        except json.JSONDecodeError as e:
            print(f"Error reading JSON file for {protein}: {e}")
            continue

        # Extract PFAM predictions and filter based on domains of interest
        for result in data.get('results', []):
            sequence_id = result.get('xref', [{}])[0].get('id')
            if not sequence_id:
                print(f"Skipping result without a sequence ID for {protein}.")
                continue

            # Extract PFAM domains from the JSON result
            sequence_pfam_domains = extract_pfam_domains(result)

            # Log the PFAM domains found for the sequence
            print(f"Domains found for {sequence_id}: {sequence_pfam_domains}")

            # Check if this sequence matches the protein based on PFAM domains
            if sequence_pfam_domains == pfam_domains_in_csv:
                protein_sequence_ids[protein].add(sequence_id)
            else:
                print(
                    f"Domains for {sequence_id} do not match CSV for {protein}. Expected: {pfam_domains_in_csv}, Found: {sequence_pfam_domains}")

        # Log the number of sequences matching the PFAM domains
        print(f"{len(protein_sequence_ids[protein])} sequences found for protein {protein}")

        # Log FASTA IDs
        fasta_ids = [record.id.split()[0] for record in SeqIO.parse(fasta_file_path, 'fasta')]
        print(f"Sequence IDs in {fasta_file_path}: {fasta_ids}")

        # Read the input FASTA file and write the filtered sequences to the corresponding output FASTA file
        try:
            with open(output_fasta_path, 'w') as output_handle:
                for record in SeqIO.parse(fasta_file_path, 'fasta'):
                    fasta_id = record.id.split()[0]  # Handle potential additional descriptions
                    if fasta_id in protein_sequence_ids[protein]:
                        SeqIO.write(record, output_handle, 'fasta')
                    else:
                        print(f"Skipping {record.id} as it doesn't match any sequence ID for {protein}")
            print(f"Filtered sequences written for {protein} to {output_fasta_path}")
        except Exception as e:
            print(f"Error processing FASTA file for {protein}: {e}")

    print("Filtered sequences written to the output directory.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process PFAM domain data and filter FASTA sequences.')
    parser.add_argument('directory', type=str, help='Directory containing the JSON and FASTA files.')
    parser.add_argument('csv_filename', type=str,
                        help='Name of the CSV file containing protein names and PFAM domains.')

    args = parser.parse_args()

    main(args.directory, args.csv_filename)
