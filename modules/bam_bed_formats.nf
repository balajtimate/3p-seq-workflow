#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CREATE_BAM {

    label 'processing'
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results/", mode: 'copy', pattern: '*'

    input:
    tuple library, bed_file, bam_file

    output:
    tuple library, path('*.bam'), emit: bam
    tuple library, path('*.bed'), emit: bed

    script:
    """
    python ${projectDir}/modules/filter_bam.py -b ${bam_file} -i ${bed_file} -o ${library}.filtered.bam -ob ${library}.filtered.bed
    """
}

process BAM_TO_BED {

    label 'processing'
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results/", mode: 'copy', pattern: '*'

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path('*.bed'), emit: bed

    script:
    """
    # Bedtools bamtobed
    bedtools bamtobed -i ${bam} -bed12 > ${library}.bed12

    # Add read length to bed12
        awk '{
        # Calculate read length as the sum of block sizes (11th column)
        split(\$11, block_sizes, ",");
        read_length = 0;
        for (i = 1; i <= length(block_sizes); i++) {
            read_length += block_sizes[i];
        }
        # Print chr, start, end, name, score, strand, and read_length
        print \$1, \$2, \$3, \$4, \$5, \$6, read_length;
    }' OFS='\t' ${library}.bed12 > ${library}.bed
    """
}