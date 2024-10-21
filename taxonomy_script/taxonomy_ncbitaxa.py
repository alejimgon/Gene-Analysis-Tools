import pandas as pd
from ete3 import NCBITaxa
import logging
import argparse

# Use: python taxonomy_ncbitaxa.py path/to/taxids.tsv path/to/blastp.txt path/to/output_file.tsv
# Input Files:
# Taxonomy File: Create a tab-separated values (TSV) file named taxids.tsv. This file should contain at least
# two columns: assembly_accession and taxid. Ensure the first row is the header.
# BLASTP Results File: Create a text file named blastp.txt. Each line should contain a unique identifier in the
# format assembly_accession-unique_id.

# Function to get taxonomy lineage
def get_taxonomy_lineage(ncbi, taxid):
    try:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        rank_names = {rank: names[taxid] for taxid, rank in ranks.items()}
        rank_order = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        return [rank_names.get(rank, '') for rank in rank_order]
    except Exception as e:
        logging.error(f"Error retrieving taxonomy for taxid {taxid}: {e}")
        return ['NA'] * 7  # Return empty ranks if there's an issue


def main(assembly_accession_taxid, blastp_result, output_file_path):
    # Define error output path and log file path
    error_file_path = 'custom_refseq_error_taxonomy.tsv'
    log_file_path = 'assembly_taxonomy_log.txt'

    # Initialize NCBITaxa
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()

    # Configure logging
    logging.basicConfig(filename=log_file_path, level=logging.ERROR, format='%(asctime)s %(message)s')

    # Read the input file
    df = pd.read_csv(assembly_accession_taxid, sep='\t')

    # Read the BLAST results
    with open(blastp_result, 'r') as f:
        unique_ids = f.read().splitlines()

    # Parse the BLAST results
    parsed_data = [line.split('-') for line in unique_ids]
    parsed_df = pd.DataFrame(parsed_data, columns=['assembly_accession', 'unique_id'])
    parsed_df['full_id'] = parsed_df['assembly_accession'] + '-' + parsed_df['unique_id']

    # Initialize lists for results and errors
    taxonomy_results = []
    error_entries = []

    # Process each row to get taxonomy information
    for _, row in df.iterrows():
        taxid = row['taxid']
        assembly_accession = row['assembly_accession']
        taxonomy = get_taxonomy_lineage(ncbi, taxid)
        if taxonomy == ['NA'] * 7:
            error_entries.append([assembly_accession, taxid])
        taxonomy_results.append([assembly_accession] + taxonomy)

    # Create DataFrames for results and errors
    taxonomy_columns = ['assembly_accession', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxonomy_df = pd.DataFrame(taxonomy_results, columns=taxonomy_columns)
    error_df = pd.DataFrame(error_entries, columns=['assembly_accession', 'taxid'])

    # Merge the parsed data with the taxonomy information
    merged_df = parsed_df.merge(taxonomy_df, on='assembly_accession', how='left')

    # Reorder columns to have 'full_id' first followed by taxonomy columns
    final_columns = ['full_id'] + taxonomy_columns[1:]
    merged_df = merged_df[final_columns]

    # Write the result to TSV files
    merged_df.to_csv(output_file_path, sep='\t', index=False)
    error_df.to_csv(error_file_path, sep='\t', index=False)

    print(f"Output written to {output_file_path}")
    print(f"Errors written to {error_file_path}")
    print(f"Log written to {log_file_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Retrieve taxonomy information based on assembly accessions.')
    parser.add_argument('assembly_accession_taxid', type=str,
                        help='Path to the TSV file containing assembly accessions and tax IDs.')
    parser.add_argument('blastp_result', type=str, help='Path to the BLASTP results text file.')
    parser.add_argument('output_file_path', type=str, help='Path to the output TSV file for taxonomy results.')

    args = parser.parse_args()

    main(args.assembly_accession_taxid, args.blastp_result, args.output_file_path)