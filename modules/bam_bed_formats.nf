#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CREATE_BAM {

    label 'processing'
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results/", mode: 'copy', pattern: '*'

    input:
    tuple val(library), path(bed_file), path(bam_file)

    output:
    tuple val(library), path('*.bam'), emit: bam
    tuple val(library), path('*.bai'), emit: bai
    tuple val(library), path('*.bed'), emit: bed

    script:
    """
    python ${projectDir}/modules/filter_bam.py -b ${bam_file} -i ${bed_file} -o ${library}.filtered.bam -ob ${library}.filtered.bed
    samtools index ${library}.filtered.bam 
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

process CREATE_BIGWIG {

    label 'processing'
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results/", mode: 'copy', pattern: '*.bw'
    publishDir "${params.out_dir}/${library}_results/", mode: 'copy', pattern: '*.bedgraph'

    input:
    tuple val(library), path(bed_file)
    path(chrom_sizes)

    output:
    tuple val(library), path('*.positive.bw'), emit: pos_bigwig
    tuple val(library), path('*.negative.bw'), emit: neg_bigwig
    tuple val(library), path('*.positive.sorted.bedgraph'), emit: pos_bedgraph
    tuple val(library), path('*.negative.sorted.bedgraph'), emit: neg_bedgraph

    script:
    """
    # Split bed file by strand
    awk '\$6 == "+" {print \$1, \$2, \$3, \$5}' OFS="\t" ${bed_file} > ${library}.positive.bedgraph
    awk '\$6 == "-" {print \$1, \$2, \$3, \$5}' OFS="\t" ${bed_file} > ${library}.negative.bedgraph

    # Sort the BedGraph files
    sort -k1,1 -k2,2n ${library}.positive.bedgraph > ${library}.positive.sorted.bedgraph
    sort -k1,1 -k2,2n ${library}.negative.bedgraph > ${library}.negative.sorted.bedgraph

    # Convert BedGraph to BigWig
    bedGraphToBigWig ${library}.positive.sorted.bedgraph ${chrom_sizes} ${library}.positive.bw
    bedGraphToBigWig ${library}.negative.sorted.bedgraph ${chrom_sizes} ${library}.negative.bw
    """
}
