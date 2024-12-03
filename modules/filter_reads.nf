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

process COUNT_CS_READS {

    label 'processing'

    tag { library }

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*"

    input:
    tuple val(library), path(bed)

    output:
    tuple val(library), path('*counts.bed'), emit: counts_cs_tuple
    path('*.3pSites.bed')

    script:
    """
    # Adjust start/end position based on read strand
    awk '{if (\$6 == "+") \$2 = \$3 - 1; else if (\$6 == "-") \$3 = \$2 + 1; print}' OFS='\t' ${bed} > ${library}.3pSites.bed

    # Reorganize columns: chr, start, end, cs_id, score, strand, length
    # Count reads at cleavage sites, aggregate by 
    awk '
    {
        # Define the key for aggregation (chromosome, start, end, strand)
        key = \$1 "\t" \$2 "\t" \$3 "\t" \$6;

        # Construct the cleavage site ID
        cleavage_site_id = \$1 ":" \$2 ":" \$3 ":" \$6;

        # Aggregate lengths and IDs for each key
        lengths[key] = lengths[key] ? lengths[key] "," \$7 : \$7;
        ids[key] = ids[key] ? ids[key] "," \$4 : \$4;

        # Store the cleavage site ID
        site_id[key] = cleavage_site_id;

        # Count occurrences of each key
        count[key]++;
    }
    END {
        for (k in count) {
            # Extract chrom, start, end, and strand for printing
            split(k, fields, "\t");
            print fields[1], fields[2], fields[3], site_id[k], count[k], fields[4], lengths[k], ids[k];
        }
    }' OFS='\t' ${library}.3pSites.bed | sort > ${library}.3pSites.counts.bed
    """
}

process FILTER_IQR {

    label 'processing'
    
    tag { library }

    conda = "${HOME}/miniconda3/envs/jupyter"

    publishDir "${params.out_dir}/${library}_results", mode: 'copy', pattern: "*"

    input:
    tuple val(library), path(bed)

    output:
    tuple val(library), path('*.filtered.bed'), emit: filtered_cs_tuple
    path("*.png")

    script:
    """
    python ${projectDir}/modules/filter_iqr.py -i ${bed} -s ${library} -o ${library}.filtered.bed
    """
}