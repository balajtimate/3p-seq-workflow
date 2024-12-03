#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CREATE_BAM {

    label 'processing'
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results/", mode: 'copy', pattern: '*'

    input:
    tuple val(library), path(bed)
    path non_overlapping_genes_bed

    output:
    tuple val(library), path('*.final.bed'), emit: unique_stranded_filtered_tuple
    // path('*.unique_mappers.bed')
    // path('*.unique_mappers.overlap.bed')
    // path('*.unique_mappers.overlap.stranded.bed')
    // path('*.unique_mappers.overlap.stranded.filtered.bed')

    script:
    """
    # Filter for unique mappers
    awk '\$5 == 255' ${bed} > ${library}.unique_mappers.bed

    # Intersect with non-overlapping genes
    bedtools intersect -a ${library}.unique_mappers.bed -b ${non_overlapping_genes_bed} -wa -wb > ${library}.unique_mappers.overlap.bed

    # Check for gene strandedness
    awk '{print \$1, \$2, \$3, \$4, \$5, \$6, \$13, \$7}' OFS='\t' ${library}.unique_mappers.overlap.bed > ${library}.unique_mappers.overlap.stranded.bed

    # Filter for stranded reads with same strand as gene
    awk '\$6 == \$7' OFS='\t' ${library}.unique_mappers.overlap.stranded.bed > ${library}.unique_mappers.overlap.stranded.filtered.bed

    # Remove reads where length >45
    awk '\$8 <= 45 {print \$1, \$2, \$3, \$4, \$5, \$7, \$8}' OFS='\t' ${library}.unique_mappers.overlap.stranded.filtered.bed > ${library}.unique_mappers.overlap.stranded.filtered.final.bed
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