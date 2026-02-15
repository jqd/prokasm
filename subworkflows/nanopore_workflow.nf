// EPI2ME-style Nanopore Assembly Workflow
nextflow.enable.dsl = 2

include { NANOPLOT }   from '../modules/qc/nanoplot'
include { BAM2FQ }     from '../utils/bam2fq.nf'
include { FLYE }       from '../modules/assembly/flye'
include { CIRCLATOR }  from '../modules/assembly/circlator'
include { DNAAPLER }   from '../modules/assembly/dnaapler'
include { BUSCO }      from '../modules/assembly/busco'
// include { MINI_ALIGN } from '../modules/report/mini_align'


/**
 * Nanopore Long-Read Assembly Workflow
 * 
 * Supports both BAM and FASTQ input formats via --longread_format parameter
 * 
 * Steps:
 * 1. BAM to FASTQ conversion (if BAM input)
 * 2. QC (NanoPlot)
 * 3. Assembly (Flye)
 * 4. Contig reorientation (dnaapler/circlator)
 * 5. Quality assessment (BUSCO)
 * 6. Read alignment for QC (minimap2)
 */
workflow NanoporeWorkflow {

    take:
    samples_ch  // tuple val(sample_id), path(reads) - BAM or FASTQ

    main:
    
    // ---------------------------------
    // Step 1: BAM to FASTQ conversion (if needed)
    // ---------------------------------
    def format = params.longread_format?.toLowerCase() ?: 'fastq'
    
    if (format == 'bam') {
        BAM2FQ(samples_ch)
        fastq_ch = BAM2FQ.out.fastq
    } else {
        // Already in FASTQ format
        fastq_ch = samples_ch
    }
    
    // ---------------------------------
    // Step 2: QC (optional)
    // ---------------------------------
    ch_qc_files = Channel.empty()
    
    if (!params.skip_qc) {
        NANOPLOT(fastq_ch)
        ch_qc_files = ch_qc_files.mix(NANOPLOT.out.stats.map { sample_id, stats -> stats })
    }

    // ---------------------------------
    // Step 3: Assembly (Flye)
    // ---------------------------------
    assembly_ch = FLYE(fastq_ch)

    // ---------------------------------
    // Step 4: Contig reorientation (optional)
    // ---------------------------------
    if (!params.skip_reorientation) {
        def reorientation_tool = params.reorientation_tool?.toLowerCase() ?: 'dnaapler'
        
        if (reorientation_tool == 'circlator') {
            reoriented_ch = CIRCLATOR(assembly_ch.contigs)
            final_contigs = reoriented_ch.contigs
        } else {
            // Default: dnaapler  
            reoriented_ch = DNAAPLER(assembly_ch.contigs)
            final_contigs = reoriented_ch.contigs
        }
    } else {
        final_contigs = assembly_ch.contigs
    }

    // ---------------------------------
    // Step 5: BUSCO (optional)
    // ---------------------------------
    if (!params.skip_busco) {
        busco_ch = BUSCO(final_contigs)
        ch_qc_files = ch_qc_files.mix(busco_ch.summary.map { sample_id, summary -> summary })
    }

    // ---------------------------------
    // Step 6: Read alignment for QC
    // Join reads with final contigs for alignment
    // ---------------------------------
    joined_ch = fastq_ch
        .join(final_contigs, by: 0)
    
    // MINI_ALIGN(joined_ch)
    // ch_qc_files = ch_qc_files.mix(MINI_ALIGN.out.stats)

    emit:
    contigs   = final_contigs
    assembly  = assembly_ch.contigs
    // bam       = MINI_ALIGN.out.bam
    // stats     = MINI_ALIGN.out.stats
    qc_files  = ch_qc_files
}
