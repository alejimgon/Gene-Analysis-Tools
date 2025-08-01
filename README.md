# Gene-Analysis-Tools

This repository contains a collection of scripts and tools for analyzing genomic data, including:

## InterPro PFAM filtering
PFAM Domain Filtering: Tools for extracting and filtering PFAM domain data from JSON and FASTA files.  
Path: pfam_filtering/interproscan_pfam_filtering.py

## Phylogenetic pipeline with Nextflow
Nextflow pipeline to compute a phylogenetic tree from FASTA files. The pipeline will align (MAFFT), trim (trimAl), and compute the phylogenetic tree (IQ-TREE) for all files in the "fasta_inputs" folder.  
Path: phylogenetic_pipeline_Nextflow/

## Presence/absence analysis and heat map 
Gene Presence/Absence Analysis: Scripts for generating heatmaps based on gene presence/absence data.  
Path: presence_absence_heat_map/heat_map.py  
Path: presence_absence_heat_map/presence_absence.py

## Taxon scaling analysis for metatranscriptomics data
Taxonomic Scaling Analysis: An R script for normalizing metatranscriptomics data for differential expression analyses.  
Path: taxon_scaling/taxon_scaling_analysis.r

## Taxonomy information
Taxonomy Retrieval: A Python script for retrieving taxonomy information using NCBI Taxonomy.  
Path: taxonomy_script/taxonomy_ncbitaxa.py

## License

This project is for non-commercial use.  
See [LICENSE](LICENSE) for details.