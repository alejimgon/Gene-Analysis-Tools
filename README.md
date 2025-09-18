# Gene-Analysis-Tools

This repository contains a collection of scripts for genomic data analysis.

## Contents

- **extract_genbank**
  Extract a region from GenBank files based on a table of gene pairs and accession numbers.  
  For each gene pair, it finds the two reference genes in the GenBank file, extracts all features between them (inclusive), plus two genes upstream and two downstream (if available), and writes the region to a new GenBank file. See [extract_genbank/README.md](extract_genbank/README.md) for usage and details.


- **extract_seq/**  
  Perl scripts to extract sequences based on a list from a FASTA file. It contains two scripts: one to extract the sequences present in a reference list and another to extract the sequences not present in a reference list. See [extract_seq/README.md](extract_seq/README.md) for usage and details.

- **pfam_filtering/**  
  PFAM domain filtering tools. See [pfam_filtering/README.md](pfam_filtering/README.md) for usage and details.

- **presence_absence_heat_map/**  
  Gene presence/absence analysis and heatmap scripts. See [presence_absence_heat_map/README.md](presence_absence_heat_map/README.md).

- **taxonomy_script/**  
  Taxonomy retrieval scripts using NCBI Taxonomy. See [taxonomy_script/README.md](taxonomy_script/README.md).

## License

This project is for non-commercial use.  
See [LICENSE](LICENSE) for details.