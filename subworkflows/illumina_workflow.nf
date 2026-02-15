nextflow.enable.dsl = 2

include { FASTQC }       from '../modules/qc/fastqc'
include { TRIMMOMATIC }  from '../modules/trim/trimmomatic'
include { FASTP }        from '../modules/trim/fastp'
include { TRIMGALORE }   from '../modules/trim/trimgalore'
include { MEGAHIT }      from '../modules/assembly/megahit'
include { METASPADES }   from '../modules/assembly/metaspades'
include { CIRCLATOR  }   from '../modules/assembly/circlator'
include { DNAAPLER }     from '../modules/assembly/dnaapler'
include { BUSCO      }   from '../modules/assembly/busco'


/**
 * Illumina Short-Read Assembly Workflow
 * 
 * Steps:
 * 1. QC (FastQC)
 * 2. Trimming (fastp/trimmomatic/trimgalore)
 * 3. Assembly (MEGAHIT/MetaSPAdes)
 * 4. Contig reorientation (dnaapler/circlator)
 * 5. Quality assessment (BUSCO)
 */
workflow IlluminaWorkflow {

    take:
    samples_ch  // tuple val(sample_id), path(R1), path(R2)

    main:
    
    // Initialize channels for tracking outputs
    assembly_ch = Channel.empty()
    reoriented_ch = Channel.empty()
    busco_ch = Channel.empty()
    
    // ---------------------------------
    // Step 1: QC (optional)
    // ---------------------------------
    ch_qc_files = Channel.empty()
    
    if (!params.skip_qc) {
        FASTQC(samples_ch)
        ch_qc_files = ch_qc_files.mix(FASTQC.out.zip.map { sample_id, zip -> zip })
    }
    
    // ---------------------------------
    // Step 2: Trimming (optional)
    // ---------------------------------
    if (!params.skip_trimming) {
        def trimmer = params.illumina_trimmer.toLowerCase()
        
        if (trimmer == 'trimmomatic') {
            trimmed_ch = TRIMMOMATIC(samples_ch)
        } else if (trimmer == 'trimgalore') {
            trimmed_ch = TRIMGALORE(samples_ch)
        } else {
            // Default: fastp
            trimmed_ch = FASTP(samples_ch)
        }
        
        reads_for_assembly = trimmed_ch.reads
    } else {
        // Skip trimming - use raw reads
        reads_for_assembly = samples_ch
    }

    // ---------------------------------
    // Step 3: Assembly
    // ---------------------------------
    def assembler = params.illumina_assembler.toLowerCase()

    if (assembler == 'metaspades') {
        assembly_ch = METASPADES(reads_for_assembly)
    } else {
        // Default: MEGAHIT
        assembly_ch = MEGAHIT(reads_for_assembly)
    }
   
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

    emit:
    contigs     = final_contigs
    assembly    = assembly_ch.contigs
    qc_files    = ch_qc_files
}
