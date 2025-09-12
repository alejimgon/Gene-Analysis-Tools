# Gene Presence/Absence Analysis & Heatmap

Scripts for generating gene presence/absence tables and heatmaps.

## Usage

### 1. Create Presence/Absence Table

```bash
python presence_absence.py path/to/fasta_directory path/to/gene_order.txt path/to/output.csv
```

### 2. Generate Heatmap

```bash
python heat_map.py path/to/output.csv path/to/species_order.txt path/to/output.pdf
```

Or use the combined script:

```bash
python pres_abs_heat_map.py --fasta_dir path/to/fasta_directory --gene_order_file path/to/gene_order.txt --species_file path/to/species_order.txt --output_csv path/to/output.csv --output_pdf path/to/output.pdf
```

## Output

- CSV file with gene presence/absence
- PDF heatmap