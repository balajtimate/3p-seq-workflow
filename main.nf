/* 
 * enables modules 
 */
nextflow.enable.dsl = 2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

log.info """\
 3'-Seq analysis - N F   P I P E L I N E
 ===================================
 """

// import modules
include { ASSIGN_STRANDEDNESS } from './modules/motif_analysis.nf'
include { COUNT_CS_READS } from './modules/motif_analysis.nf'
include { MOTIF_ANALYSIS } from './modules/motif_analysis.nf'

genome_fa_ch = Channel.fromPath(params.genome_fa, checkIfExists: true).collect()
non_overlap_genes_ch = Channel.fromPath(params.non_overlap_genes, checkIfExists: true).collect()
polya_sites_bed_ch = Channel.fromPath(params.polya_sites_bed, checkIfExists: true).collect()

// Subworkflow for preprocessing steps
workflow motif_analysis {
    take:
        input_bed

    main:
        ASSIGN_STRANDEDNESS(input_bed, non_overlap_genes_ch)
        unique_stranded_filtered_tuple = ASSIGN_STRANDEDNESS.out.unique_stranded_filtered_tuple
        MOTIF_ANALYSIS(unique_stranded_filtered_tuple, genome_fa_ch, polya_sites_bed_ch)
        motif_tuple = MOTIF_ANALYSIS.out.motif_tuple
        // COUNT_CS_READS(unique_stranded_filtered_tuple)
        // counts_cs_tuple = COUNT_CS_READS.out.counts_cs_tuple
    emit:
        motif_tuple
}

// Main workflow
workflow {
    if (params.run_mode == 'motif_analysis') {
        input_bed_ch = Channel.fromPath(params.input_bed, checkIfExists: true).map { input_bed_path -> tuple(input_bed_path.baseName, input_bed_path) }
        input_bed_ch.each {
            motif_analysis(input_bed_ch)
        }
    }
}

/* 
 * completion handler
 */
workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}