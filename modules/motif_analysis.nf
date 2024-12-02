#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process ASSIGN_STRANDEDNESS {

    label 'processing'
    
    tag { library }

    publishDir "${params.out_dir}/${library}_results/", mode: 'copy', pattern: '*'

    input:
    tuple val(library), path(bed)
    path non_overlapping_genes_bed

    output:
    tuple val(library), path('*.final.bed'), emit: unique_stranded_filtered_tuple

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

process COUNT_CS_READS {

    label 'processing'

    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*.bed"

    input:
    tuple val(library), path(bed)

    output:
    tuple val(library), path('*counts.bed'), emit: counts_cs_tuple

    script:
    """
    # Count reads at cleavage sites, aggregate by 
    awk '
    {
        key = \$1 "\t" \$2 "\t" \$3 "\t" \$6;
        lengths[key] = lengths[key] ? lengths[key] "," \$7 : \$7;
        ids[key] = ids[key] ? ids[key] "," \$4 : \$4;
        count[key]++;
    }
    END {
        for (k in count) {
            print k, count[k], lengths[k], ids[k];
        }
    }' OFS='\t' ${bed} | sort > ${library}.3pSites.counts.bed
    """
}

process MOTIF_ANALYSIS {

    label 'motif_analysis'

    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*canonical.bed"
    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*.3pSites.bed"  

    input:
    tuple val(library), path(bed)
    path genome_fa
    path canon_polya_sites_bed

    output:
    tuple val(library), path('*canonical.bed'), emit: motif_tuple
    path('*.3pSites.bed')

    script:
    """
    # Adjust start/end position based on read strand
    awk '{if (\$6 == "+") \$2 = \$3 - 1; else if (\$6 == "-") \$3 = \$2 + 1; print}' OFS='\t' ${bed} > ${library}.3pSites.bed

    # Count AAUAAA motif
    python ${projectDir}/modules/check_AAUAAA_motif.py -i ${library}.3pSites.bed -f ${genome_fa} -o ${library}.AAUAAA.bed

    # Count canonical sites
    python ${projectDir}/modules/count_canon_sites.py -i ${library}.AAUAAA.bed -c ${canon_polya_sites_bed} -o ${library}.AAUAAA.canonical.bed
    """
}
