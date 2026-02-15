nextflow.enable.dsl = 2

include { BAM2FQ }           from '../utils/bam2fq.nf'
include { UNICYCLER_HYBRID } from '../modules/assembly/unicycler_hybrid.nf'
include { FASTQC }           from '../modules/qc/fastqc.nf'
include { NANOPLOT }         from '../modules/qc/nanoplot.nf'
include { FASTP }            from '../modules/trim/fastp.nf'
include { TRIMMOMATIC }      from '../modules/trim/trimmomatic'
include { TRIMGALORE }       from '../modules/trim/trimgalore'
include { CIRCLATOR }        from '../modules/assembly/circlator'
include { DNAAPLER }         from '../modules/assembly/dnaapler'
include { BUSCO }            from '../modules/assembly/busco'


/**
 * Hybrid Assembly Workflow
 * 
 * Combines pre-merged Illumina short reads with long reads (PacBio or Nanopore)
 * for hybrid genome assembly using Unicycler
 * 
 * Steps:
 * 1. BAM to FASTQ conversion (if long reads are in BAM format)
 * 2. Illumina QC (FastQC)
 * 3. Illumina trimming (fastp/trimmomatic/trimgalore)
 * 4. Long read QC (NanoPlot)
 * 5. Hybrid assembly (Unicycler)
 * 6. Contig reorientation (dnaapler/circlator)
 * 7. Quality assessment (BUSCO)
 * 
 * Input:
 *   samplesheet: Channel of [sample_id, R1, R2, longreads] tuples (already merged)
 *   longread_type: String - platform type (pacbio or nanopore)
 *   longread_format: String - format type (bam or fastq)
 * 
 * Output:
 *   contigs: Assembled and reoriented contigs from Unicycler
 */
workflow HybridWorkflow {
    take: 
        samplesheet       // [sample_id, R1, R2, longreads]
        longread_type     // pacbio or nanopore
        longread_format   // bam or fastq
    
    main:
        log.info "Processing Hybrid Assembly with merged reads..."
        log.info "Long read type: ${longread_type}, format: ${longread_format}"
        
        // ---------------------------------
        // Step 1: Convert BAM to FASTQ if needed
        // ---------------------------------
        if (longread_format == 'bam') {
            log.info "Converting BAM longreads to FASTQ..."
            
            longread_ch = samplesheet.map { sample_id, illumina_R1, illumina_R2, longreads ->
                [sample_id, longreads]
            }
            BAM2FQ(longread_ch)
            
            // Recombine with illumina reads
            illumina_ch = samplesheet.map { sample_id, illumina_R1, illumina_R2, longreads ->
                [sample_id, illumina_R1, illumina_R2]
            }
            
            hybrid_data = illumina_ch.join(BAM2FQ.out.fastq)
        } else {
            // Longreads already in FASTQ format
            hybrid_data = samplesheet
        }
        
        // Split channels for parallel QC/trimming
        illumina_ch = hybrid_data.map { sample_id, illumina_R1, illumina_R2, longreads ->
            [sample_id, illumina_R1, illumina_R2]
        }
        
        longread_ch = hybrid_data.map { sample_id, illumina_R1, illumina_R2, longreads ->
            [sample_id, longreads]
        }

        // ---------------------------------
        // Step 2: Illumina QC (optional)
        // ---------------------------------
        ch_qc_files = Channel.empty()
        
        if (!params.skip_qc) {
            FASTQC(illumina_ch)
            ch_qc_files = ch_qc_files.mix(FASTQC.out.zip.map { sample_id, zip -> zip })
        }
        
        // ---------------------------------
        // Step 3: Illumina trimming (optional)
        // ---------------------------------
        if (!params.skip_trimming) {
            def trimmer = params.illumina_trimmer.toLowerCase()

            if (trimmer == 'trimmomatic') {
                trimmed_ch = TRIMMOMATIC(illumina_ch)
            } else if (trimmer == 'trimgalore') {
                trimmed_ch = TRIMGALORE(illumina_ch)
            } else {
                // Default: fastp
                trimmed_ch = FASTP(illumina_ch)
            }
            
            trimmed_reads = trimmed_ch.reads
        } else {
            // Skip trimming - use raw reads
            trimmed_reads = illumina_ch
        }
        
        // ---------------------------------
        // Step 4: Long read QC (optional)
        // ---------------------------------
        if (!params.skip_qc) {
            NANOPLOT(longread_ch)
            ch_qc_files = ch_qc_files.mix(NANOPLOT.out.stats.map { sample_id, stats -> stats })
        }
        
        // ---------------------------------
        // Step 5: Unicycler hybrid assembly
        // ---------------------------------
        // Recombine trimmed Illumina with longreads for assembly
        assembly_input = trimmed_reads.join(longread_ch)
        
        UNICYCLER_HYBRID(assembly_input)
        assembly_ch = UNICYCLER_HYBRID.out

        // ---------------------------------
        // Step 6: Contig reorientation (optional)
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
        // Step 7: BUSCO (optional)
        // ---------------------------------
        if (!params.skip_busco) {
            busco_ch = BUSCO(final_contigs)
            ch_qc_files = ch_qc_files.mix(busco_ch.summary.map { sample_id, summary -> summary })
        }
    
    emit:
    contigs  = final_contigs
    assembly = assembly_ch.contigs
    qc_files = ch_qc_files
}