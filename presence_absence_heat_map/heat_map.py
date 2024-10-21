import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# Use python gene_presence_absence_heatmap.py path/to/input.csv path/to/species_order.txt path/to/output.pdf
# Inputs:
# CSV File: a CSV file containing gene presence/absence data. The first column should contain species names
# (row indices), and subsequent columns should contain gene presence (1) or absence (0) indicators.
# Use presence_absence.py to create this file
# Species Order File: Create a text file that specifies the desired order of species, with one species name per line.

def load_species_order(species_file):
    """Load species order from a file."""
    with open(species_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def main(input_csv, species_file, output_pdf):
    # Load the presence/absence table
    df = pd.read_csv(input_csv, index_col=0)

    # Load the desired order of species (rows)
    species_order = load_species_order(species_file)

    # Ensure all specified species are in the dataframe
    for species in species_order:
        if species not in df.index:
            df.loc[species] = 0

    # Reorder the dataframe according to the specified order
    df = df.reindex(species_order)

    # Plot the heatmap
    plt.figure(figsize=(10, 15))  # Adjust figure size as needed
    heatmap = sns.heatmap(df, cmap="YlGnBu", linewidths=0.5, linecolor='gray', cbar_kws={"shrink": .5})

    # Customize the heatmap
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, fontsize=8)  # Keep species names horizontal
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=90,
                            fontsize=8)  # Rotate gene names for better readability

    # Save the heatmap as a PDF
    plt.tight_layout()
    plt.savefig(output_pdf, format='pdf')
    print(f"Heatmap saved as '{output_pdf}'.")

    # Show the plot
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a heatmap of gene presence/absence data.')
    parser.add_argument('input_csv', type=str, help='Input CSV file containing gene presence/absence data.')
    parser.add_argument('species_file', type=str,
                        help='File containing the desired species order (one species per line).')
    parser.add_argument('output_pdf', type=str, help='Output path for the heatmap PDF file.')

    args = parser.parse_args()

    main(args.input_csv, args.species_file, args.output_pdf)