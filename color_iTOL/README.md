# iTOL Taxonomy Coloring Scripts

Scripts for automatically generating iTOL (Interactive Tree of Life) dataset files with taxonomic coloring and labels for phylogenetic trees.

## Scripts

- **generate_itol_taxonomy_colors.py**  
  Generate iTOL DATASET_RANGE and DATASET_TEXT files for a single phylogenetic tree with automatic taxonomic coloring based on NCBI taxonomy.

- **batch_generate_itol.py**  
  Batch process multiple tree files in a directory to generate iTOL files for all trees at once.

## Usage

### Single Tree

```bash
python generate_itol_taxonomy_colors.py <tree_file> <taxonomy_file> -o <output_file> -l phylum --labels --save-rooted -v
```

**Arguments:**
- `tree_file`: Path to Newick tree file
- `taxonomy_file`: Path to taxonomy TSV file
- `-o, --output`: Output iTOL file path (default: auto-generated)
- `-l, --level`: Taxonomic level for grouping (default: phylum)
- `--labels`: Generate external text labels file
- `--save-rooted`: Save the rooted tree to a new file
- `-v, --verbose`: Verbose output
- `--root-method`: Tree rooting method (midpoint, outgroup, none)
- `--outgroup`: Outgroup sequence name (required for outgroup rooting)

### Batch Processing

```bash
python batch_generate_itol.py
```

**Before running, edit the script to set:**
- `trees_dir`: Directory containing `.treefile` files
- `taxonomy_file`: Path to taxonomy TSV file
- `output_dir`: Output directory for generated files

## Output Files

For each tree, the scripts generate:
- `{gene}_itol_colors.txt`: iTOL DATASET_RANGE file with taxonomic coloring
- `{gene}_itol_labels.txt`: iTOL DATASET_TEXT file with external labels
- `{gene}_midpoint.treefile`: Rooted tree file (if `--save-rooted` is used)

## Input File Format

**Taxonomy TSV:** Must include columns:
- `full_id`: Sequence identifier matching tree leaf names
- `superkingdom`, `phylum`, `class`, `order`, `family`: Taxonomic ranks

## Requirements

- Python ≥3.7
- pandas
- biopython

Install with:
```bash
pip install -r requirements.txt
```

## Notes

- The script identifies monophyletic clusters within taxonomic groups
- Colors are assigned based on phylum abundance in my personal databases (change if needed)
- Rare phyla are grouped as "Other Phylum" in the legend
- Default rooting method is midpoint rooting for consistent topology