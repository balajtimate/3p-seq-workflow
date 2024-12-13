#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process ANALYZE_MOTIFS {

    label 'motif_analysis'

    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*.bed"

    input:
    tuple val(library), path(bed)
    path genome_fa

    output:
    tuple val(library), path('*.motifs.bed'), emit: motif_tuple
    tuple val(library), path('*.aaaa.bed')

    script:
    """
    # Analyze AAAA stretches
    python ${projectDir}/modules/check_AAAA.py -i ${bed} -f ${genome_fa} -o ${library}.aaaa.bed

    # Use the output of AAAA stretch analysis to analyze motifs
    python ${projectDir}/modules/check_motifs.py -i ${library}.aaaa.bed -f ${genome_fa} -o ${library}.motifs.bed
    """
}