#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Align all FASTA files in 'fasta_inputs' using MAFFT.
 * Output files are named <original>.mafft and copied to 'results/mafft'.
 */

process ALIGNMENT {
    publishDir 'results/mafft', mode: 'copy'

    input:
        file fasta_file

    output:
        file "${fasta_file.simpleName}.mafft"

    script:
    """
    mafft ${fasta_file} > ${fasta_file.simpleName}.mafft
    """
}

/*
 * Trimming all MAFFT output files using trimAl.
 * Output files are named <original>.trimal and copied to 'results/trimal'.
 */

process TRIMMING {
    publishDir 'results/trimal', mode: 'copy'

    input:
        file mafft_file

    output:
        file "${mafft_file.simpleName}.trimal"

    script:
    """
    trimal -in ${mafft_file} -out ${mafft_file.simpleName}.trimal -gt 0.5
    """
}

process PHYLO_TREE {
    publishDir 'results/iqtree', mode: 'copy'

    input:
        file trimal_file

    output:
        file "${trimal_file}.treefile"

    script:
    """
    iqtree2 -s ${trimal_file} -m MFP -bb 1000 --alrt 1000 -T AUTO
    """

} 

// Define the main workflow
workflow {
    // Create a channel with all files in 'fasta_inputs'
    fasta_ch = Channel.fromPath('fasta_inputs/*')

    // Run MAFFT alignment for each file
    aligned = ALIGNMENT(fasta_ch)
    trimmed = TRIMMING(aligned)
    PHYLO_TREE(trimmed)
  
}