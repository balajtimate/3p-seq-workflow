#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process ANALYZE_MOTIFS {

    label 'motif_analysis'

    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*.motifs.bed"

    input:
    tuple val(library), path(bed)
    path genome_fa

    output:
    tuple val(library), path('*.motifs.bed'), emit: motif_tuple

    script:
    """
    # Count AAUAAA motif
    python ${projectDir}/modules/check_motifs.py -i ${bed} -f ${genome_fa} -o ${library}.motifs.bed
    """
}