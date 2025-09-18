# Extract Short GenBank Regions

This script extracts a region from GenBank files based on a table of gene pairs and accession numbers.  
For each gene pair, it finds the two reference genes in the GenBank file, extracts all features between them (inclusive), plus two genes upstream and two downstream (if available), and writes the region to a new GenBank file.

## Usage

1. **Prepare your input files:**
   - A TSV file with three columns per line:  
     ```
     gene1    gene2    accession
     ```
     Example:
     ```
     geneA    geneB    NC_000000
     ```

   - A directory containing GenBank files named as `<accession>.gbff` (e.g., `NC_000000.gbff`).

2. **Edit the script paths if needed:**
   - Set the correct paths for your TSV file, GenBank directory, and output directory at the top of the script:
     ```python
     TSV_PATH = BASE_DIR / "path/to/gene_pairs.tsv"
     GENBANK_DIR = BASE_DIR / "path/to/genbank_files"
     OUTPUT_DIR = BASE_DIR / "path/to/short_genbanks"
     ```

3. **Run the script:**
   ```bash
   python extract_short_genbank.py
   ```

## Output

- For each line in the TSV, a new GenBank file will be created in the output directory, named:
  ```
  <accession>_short_<gene1>_<gene2>.gbk
  ```
- Each output file contains the region between the two genes (inclusive), plus two genes upstream and downstream if available.

## Requirements

- Python 3
- [Biopython](https://biopython.org/) (`pip install biopython`)

## Notes

- If either gene is not found in the GenBank file, a warning will be printed and that entry will be skipped.
- The script assumes GenBank files are named as `<accession>.gbff` and located in the specified directory.

---
```# Extract Short GenBank Regions

This script extracts a region from GenBank files based on a table of gene pairs and accession numbers.  
For each gene pair, it finds the two reference genes in the GenBank file, extracts all features between them (inclusive), plus two genes upstream and two downstream (if available), and writes the region to a new GenBank file.

## Usage

1. **Prepare your input files:**
   - A TSV file with three columns per line:  
     ```
     gene1    gene2    accession
     ```
     Example:
     ```
     geneA    geneB    NC_000000
     ```

   - A directory containing GenBank files named as `<accession>.gbff` (e.g., `NC_000000.gbff`).

2. **Edit the script paths if needed:**
   - Set the correct paths for your TSV file, GenBank directory, and output directory at the top of the script:
     ```python
     TSV_PATH = BASE_DIR / "path/to/gene_pairs.tsv"
     GENBANK_DIR = BASE_DIR / "path/to/genbank_files"
     OUTPUT_DIR = BASE_DIR / "path/to/short_genbanks"
     ```

3. **Run the script:**
   ```bash
   python extract_short_genbank.py
   ```

## Output

- For each line in the TSV, a new GenBank file will be created in the output directory, named:
  ```
  <accession>_short_<gene1>_<gene2>.gbk
  ```
- Each output file contains the region between the two genes (inclusive), plus two genes upstream and downstream if available.

## Requirements

- Python 3
- [Biopython](https://biopython.org/) (`pip install biopython`)

## Notes

- If either gene is not found in the GenBank file, a warning will be printed and that entry will be skipped.
- The script assumes GenBank files are named as `<accession>.gbff` and located in the specified directory.

---