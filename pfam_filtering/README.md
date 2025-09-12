# PFAM Filtering Tools

Scripts for extracting and filtering PFAM domain data from InterProScan JSON and FASTA files.

## Usage

```bash
python interproscan_pfam_filtering.py /path/to/folder protein_pfam.csv
```

- **JSON Files:** Each protein should have a `<protein_name>.json` file.
- **FASTA Files:** Each protein should have a `<protein_name>.faa` file.
- **CSV File:** Contains protein names and their PFAM domains.

## Output

Filtered FASTA files for each protein, containing only sequences matching the specified PFAM domains.