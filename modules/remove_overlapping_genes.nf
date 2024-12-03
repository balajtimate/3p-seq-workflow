#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process REMOVE_OVERLAPPING_GENES {

    executor = 'slurm'
    clusterOptions = "--qos=30min"
    time = { 30.m }
    memory = { 8.GB }
    cpus = 8
    conda = "${HOME}/miniconda3/envs/ipa-immune"

    publishDir "${params.out_dir}", mode: 'copy', pattern: '*'

    input:
    path(gtf)

    output:
    path('*.non_overlapping_genes.bed'), emit: non_overlapping_genes_bed

    script:
    """
    filename=\$(basename ${gtf})
    gtfname="\${filename%.*}"

    awk '\$3 == "gene" && \$0 ~ /gene_type "protein_coding"/' ${gtf} > \$gtfname.pcg.gtf

    gtf2bed < \$gtfname.pcg.gtf > \$gtfname.pcg.bed

    bedtools sort -i \$gtfname.pcg.bed > \$gtfname.pcg.sorted.bed

    bedtools merge -i \$gtfname.pcg.sorted.bed -c 4 -o count | awk '\$4 == 1' | bedtools intersect -a \$gtfname.pcg.sorted.bed -b - -wa > \$gtfname.non_overlapping_genes.bed
    """
}

workflow {
    REMOVE_OVERLAPPING_GENES(params.gtf)
}