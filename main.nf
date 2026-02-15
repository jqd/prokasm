nextflow.enable.dsl=2

include { validateParameters; paramsSummaryLog; paramsHelp; samplesheetToList } from 'plugin/nf-schema'

include { IlluminaWorkflow }  from './subworkflows/illumina_workflow.nf'
include { NanoporeWorkflow }  from './subworkflows/nanopore_workflow.nf'
include { PacBioWorkflow }    from './subworkflows/pacbio_workflow.nf'
include { HybridWorkflow }    from './subworkflows/hybrid_workflow.nf'
include { CAT_FASTQ }         from './utils/CAT_FASTQ.nf'
include { CAT_LONGREAD }      from './utils/CAT_LONGREAD.nf'
include { MULTIQC }           from './modules/qc/multiqc.nf'
include { GET_SOFTWARE_VERSIONS } from './modules/report/software_versions.nf'
include { banner }            from './utils/help.nf'


//
// unknown recognition error type: groovyjarjarantlr4.v4.runtime.LexerNoViableAltException
// A common cause is having a space character immediately after a backslash (\) used for line continuation in a Bash script.
// Solution: Ensure there are no spaces after the backslash character when continuing a line.
//

// Reference links for development:
// https://github.com/nf-core/bacass
// https://github.com/rrwick/Perfect-bacterial-genome-tutorial
// https://github.com/EBI-Metagenomics/mettannotator
// https://github.com/epi2me-labs/wf-bacterial-genomes/tree/master
// https://github.com/galaxyproject/tools-iuc/tree/master/tools/unicycler

def helpMessage() {
    banner()
    log.info """
    DESCRIPTION:
    Nextflow DSL2 pipeline for de novo genome assembly from Illumina, PacBio, Nanopore, and hybrid data.
    Supports short-read, long-read, and hybrid approaches (primarily Illumina + Nanopore).
    Includes QC, trimming, assembly, circularization, and quality assessment (BUSCO).
    
    Supports multiple lanes/files per sample for all data types.
    
    USAGE:
    nextflow run main.nf --data_type <TYPE> --samplesheet <FILE> --outdir <path> [OPTIONS]
    
    MANDATORY ARGUMENTS:
        --outdir <path>
            Output directory (default: ./results)

        --data_type <illumina|pacbio|nanopore|hybrid>
            Type of sequencing data to process
            illumina: Short read Illumina data
            pacbio:   PacBio long read data
            nanopore: Nanopore long read data
            hybrid:   Combined short and long reads (requires --longread_type and --longread_format)
            
        --samplesheet <file>
            CSV file with sample information (columns: sample_id,R1,R2,LReads)
            Format by data type:
            illumina: Fill R1 and R2, leave LReads empty
            pacbio:   Fill LReads (or FASTQ file), leave R1 and R2 empty
            nanopore: Fill LReads (or FASTQ file), leave R1 and R2 empty
            hybrid:   Fill R1, R2, and LReads (longreads as BAM or FASTQ)
            
            Multi-lane support: List same sample_id multiple times with different files

    CONDITIONAL ARGUMENTS (data_type specific):
        --longread_type <pacbio|nanopore>
            Required when --data_type hybrid
            Specifies long read platform for hybrid assemblies
        
        --longread_format <bam|fastq>
            Required for PacBio, Nanopore, and Hybrid workflows.
            Specifies long read file format.
            bam:   Long reads in BAM format (will be converted to FASTQ)
            fastq: Long reads in FASTQ/FASTQ.GZ format
        
        --flye_mode <mode>
            Required when --data_type pacbio, nanopore
            Assembly mode when using Flye assembler
            Options: nano-raw, nano-corr, nano-hq, pacbio-raw, pacbio-corr, pacbio-hifi
    
    TOOL SELECTION:
        --illumina_trimmer <fastp|trimmomatic|trimgalore>
            Trimming tool for Illumina data (default: fastp)
            
        --illumina_assembler <megahit|metaspades>
            Assembly tool for Illumina data (default: megahit)
            
        --pacbio_assembler <hifiasm|flye>
            Assembly tool for PacBio data (default: hifiasm)
        
        --reorientation_tool <circlator|dnaapler>
            Tool for contig reorientation (default: dnaapler)
    
    SKIP OPTIONS:
        --skip_qc                Skip QC steps (FastQC/NanoPlot)
        --skip_trimming          Skip trimming for Illumina reads
        --skip_polish            Skip polishing step
        --skip_reorientation     Skip contig reorientation
        --skip_busco             Skip BUSCO assessment
    
    TOOL-SPECIFIC OPTIONS:
        --hifiasm_k <kmer_size>
            K-mer size for Hifiasm assembler (default: 21)
        --hifiasm_a <number>
            Minimum overlap length for Hifiasm assembler (default: 10)
        --unicycler_mode <mode>
            Assembly mode for Unicycler (default: normal)
        --unicycler_min_fasta_length <length>
            Minimum contig length for Unicycler output (default: 100)
        --metaspades_k <kmer_sizes>
            K-mer sizes for MetaSPAdes assembler (default: auto)
        --circlator_min_id <percentage>
            Minimum identity percentage for Circlator circularization (default: 90)
        --trimgalore_quality <value>
            Quality threshold for TrimGalore (default: 20)
        --trimgalore_stringency <value>
            Stringency for adapter removal in TrimGalore (default: 1)
        --busco_mode <mode>
            Mode for BUSCO analysis (default: genome)
    
    OPTIONAL PARAMETERS AND ADVANCED USAGE:
    Each tool supports custom parameters via --{tool}_opts flags:
        --trimmomatic_opts ""     Example: "ILLUMINACLIP:adapters.fa:2:30:10"
        --fastp_opts ""           Example: "--cut_front --cut_tail -l 50"
        --trimgalore_opts ""      Example: "--paired --illumina"
        --fastqc_opts ""          Example: "--nogroup"
        --megahit_opts ""         Example: "--min-contig-len 500"
        --spades_opts ""          Example: "--careful -k 21,33,55"
        --nanoplot_opts ""        Example: "--plots dot --legacy"
        --hifiasm_opts ""         Example: "-l 3"
        --flye_opts ""            Example: "--genome-size 5m --asm-coverage 50"
        --circlator_opts ""       Example: "--merge_min_id 85"
        --dnaapler_opts ""        Example: "--evalue 1e-10"
        --busco_opts ""           Example: "--offline"
        --unicycler_opts ""       Example: "--no_correct"
    For a full list of options for each tool, refer to the respective tool documentation.

    EXAMPLES:
    # Illumina assembly with MEGAHIT
    nextflow run main.nf --data_type illumina --samplesheet illumina.csv -profile cluster -resume --outdir illumina_results


    # PacBio HiFi assembly (BAM input)
    nextflow run main.nf --data_type pacbio --samplesheet pacbio.csv --flye_mode pacbio-hifi --longread_format fastq -profile cluster -resume --outdir pacbio_results

    # Nanopore assembly with high-quality reads
    nextflow run main.nf --data_type nanopore --samplesheet nanopore.csv --flye_mode nano-hq --longread_format fastq -profile cluster -resume --outdir nanopore_results
    
    # Hybrid assembly (Illumina + Nanopore)
    nextflow run main.nf --data_type hybrid --samplesheet hybrid.csv --longread_type nanopore --longread_format fastq -profile cluster -resume --outdir hybrid_results

    For more information, see documentation or run with --help
    """
}

if ( params.help) {
  helpMessage()
  exit 0
}

// Validate parameters using nf-schema
validateParameters()

// Log parameter summary
log.info paramsSummaryLog(workflow)

// ----------------------------
// Helper functions
// ----------------------------

/**
 * Check if a samplesheet field has a valid value
 */
def hasValidValue(value) {
    return value && value != [] && value != ''
}

// ----------------------------
// Workflow definition
// ----------------------------

workflow {
    // Parse samplesheet using nf-schema with unified schema
    def data_type = params.data_type.toLowerCase()
    
    // Use single schema for all data types
    def schema_file = "${projectDir}/assets/schema_samplesheet.json"
    
    // Parse and validate samplesheet using nf-schema
    // Note: nf-schema will validate file existence due to "format": "file-path" in schema
    samplesheet = Channel.fromList(samplesheetToList(params.samplesheet, schema_file))
        .ifEmpty { error("Samplesheet is empty or invalid: ${params.samplesheet}") }
        .map { row ->
            // Convert list to map based on schema column order
            return [
                sample_id: row[0],
                R1: row[1],
                R2: row[2],
                LReads: row[3]
            ]
        }
        .map { row ->
            // Helper check functions
            def hasR1R2 = hasValidValue(row.R1) && hasValidValue(row.R2)
            def hasLReads = hasValidValue(row.LReads)
            
            // Validate columns match data_type
            switch (data_type) {
                case 'illumina':
                    if (!hasR1R2) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: Invalid Illumina samplesheet                      ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Issue: Missing R1 and/or R2 columns
                            
                            For --data_type illumina, both R1 and R2 must be provided.
                            
                            Expected format:
                              sample_id,R1,R2,LReads
                              sample1,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,
                            
                            Samplesheet: ${params.samplesheet}
                        """)
                    }
                    if (hasLReads) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: Invalid Illumina samplesheet                      ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Issue: LReads column should be empty for Illumina data
                            
                            For --data_type illumina, only R1 and R2 should be filled.
                            Leave the LReads column empty.
                            
                            Samplesheet: ${params.samplesheet}
                        """)
                    }
                    return [row.sample_id, file(row.R1), file(row.R2)]
                    
                case 'pacbio':
                    if (!hasLReads) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: Invalid PacBio samplesheet                        ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Issue: Missing LReads column
                            
                            For --data_type pacbio, the LReads column must be provided.
                            
                            Expected format:
                              sample_id,R1,R2,LReads
                              sample1,,,/path/to/reads.bam
                            
                            Don't forget to specify --longread_format (bam or fastq)
                            
                            Samplesheet: ${params.samplesheet}
                        """)
                    }
                    if (hasR1R2) {
                        log.warn("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  WARNING: R1/R2 ignored for PacBio assembly               ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            
                            For --data_type pacbio, only LReads will be used for assembly.
                            R1 and R2 files will be IGNORED.
                            
                            If you want to use both short and long reads together,
                            change --data_type to 'hybrid' and specify --longread_type pacbio
                            
                            Samplesheet: ${params.samplesheet}
                        """)
                    }
                    
                    // Validate file format matches longread_format parameter
                    def longread_file = row.LReads.toString().toLowerCase()
                    def format_param = params.longread_format.toLowerCase()
                    
                    if (format_param == 'bam' && !longread_file.endsWith('.bam')) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: File format mismatch (PacBio)                     ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Parameter: --longread_format bam
                            File: ${row.LReads}
                            
                            The file does not have a .bam extension.
                            
                            Either:
                              1. Change --longread_format to 'fastq', OR
                              2. Use a BAM file
                        """)
                    }
                    if (format_param == 'fastq' && !longread_file.matches('.*\\.(fq|fastq|fq\\.gz|fastq\\.gz)$')) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: File format mismatch (PacBio)                     ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Parameter: --longread_format fastq
                            File: ${row.LReads}
                            
                            The file does not have a FASTQ extension (.fq, .fastq, .fq.gz, .fastq.gz).
                            
                            Either:
                              1. Change --longread_format to 'bam', OR
                              2. Use a FASTQ file
                        """)
                    }
                    
                    return [row.sample_id, file(row.LReads)]
                    
                case 'nanopore':
                    if (!hasLReads) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: Invalid Nanopore samplesheet                      ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Issue: Missing LReads column
                            
                            For --data_type nanopore, the LReads column must be provided.
                            
                            Expected format:
                              sample_id,R1,R2,LReads
                              sample1,,,/path/to/reads.fastq.gz
                            
                            Don't forget to specify --longread_format and --flye_mode
                            
                            Samplesheet: ${params.samplesheet}
                        """)
                    }
                    if (hasR1R2) {
                        log.warn("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  WARNING: R1/R2 ignored for Nanopore assembly             ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            
                            For --data_type nanopore, only LReads will be used for assembly.
                            R1 and R2 files will be IGNORED.
                            
                            If you want to use both short and long reads together,
                            change --data_type to 'hybrid' and specify --longread_type nanopore
                            
                            Samplesheet: ${params.samplesheet}
                        """)
                    }
                    
                    // Validate file format matches longread_format parameter
                    def longread_file = row.LReads.toString().toLowerCase()
                    def format_param = params.longread_format.toLowerCase()
                    
                    if (format_param == 'bam' && !longread_file.endsWith('.bam')) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: File format mismatch (Nanopore)                   ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Parameter: --longread_format bam
                            File: ${row.LReads}
                            
                            The file does not have a .bam extension.
                            
                            Either:
                              1. Change --longread_format to 'fastq', OR
                              2. Use a BAM file
                        """)
                    }
                    if (format_param == 'fastq' && !longread_file.matches('.*\\.(fq|fastq|fq\\.gz|fastq\\.gz)$')) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: File format mismatch (Nanopore)                   ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Parameter: --longread_format fastq
                            File: ${row.LReads}
                            
                            The file does not have a FASTQ extension (.fq, .fastq, .fq.gz, .fastq.gz).
                            
                            Either:
                              1. Change --longread_format to 'bam', OR
                              2. Use a FASTQ file
                        """)
                    }
                    
                    return [row.sample_id, file(row.LReads)]
                    
                case 'hybrid':
                    if (!hasR1R2) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: Invalid Hybrid samplesheet                        ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Issue: Missing R1 and/or R2 columns
                            
                            For --data_type hybrid, both R1 and R2 must be provided.
                            
                            Expected format:
                              sample_id,R1,R2,LReads
                              sample1,/path/to/R1.fq.gz,/path/to/R2.fq.gz,/path/to/longreads.fastq.gz
                            
                            Don't forget: --longread_type and --longread_format are required!
                            
                            Samplesheet: ${params.samplesheet}
                        """)
                    }
                    if (!hasLReads) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: Invalid Hybrid samplesheet                        ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Issue: Missing LReads column
                            
                            For --data_type hybrid, LReads column must be provided.
                            
                            Expected format:
                              sample_id,R1,R2,LReads
                              sample1,/path/to/R1.fq.gz,/path/to/R2.fq.gz,/path/to/longreads.fastq.gz
                            
                            Don't forget: --longread_type and --longread_format are required!
                            
                            Samplesheet: ${params.samplesheet}
                        """)
                    }
                    
                    // Validate longread_format matches actual file extension
                    def longread_file = row.LReads.toString().toLowerCase()
                    def format_param = params.longread_format.toLowerCase()
                    
                    if (format_param == 'bam' && !longread_file.endsWith('.bam')) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: File format mismatch (Hybrid)                     ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Parameter: --longread_format bam
                            File: ${row.LReads}
                            
                            The file does not have a .bam extension.
                            
                            Either:
                              1. Change --longread_format to 'fastq', OR
                              2. Use a BAM file
                        """)
                    }
                    if (format_param == 'fastq' && !longread_file.matches('.*\\.(fq|fastq|fq\\.gz|fastq\\.gz)$')) {
                        error("""
                            ╔═══════════════════════════════════════════════════════════╗
                            ║  ERROR: File format mismatch (Hybrid)                     ║
                            ╚═══════════════════════════════════════════════════════════╝
                            
                            Sample: ${row.sample_id}
                            Parameter: --longread_format fastq
                            File: ${row.LReads}
                            
                            The file does not have a FASTQ extension (.fq, .fastq, .fq.gz, .fastq.gz).
                            
                            Either:
                              1. Change --longread_format to 'bam', OR
                              2. Use a FASTQ file
                        """)
                    }
                    
                    return [row.sample_id, file(row.R1), file(row.R2), file(row.LReads)]
                    
                default:
                    error("Unknown data_type: ${data_type}")
            }
        }
    
    // =========================================
    // Handle multi-lane merging for each data type
    // =========================================
    
    // --- ILLUMINA: Merge multiple lanes ---
    if (data_type == 'illumina') {
        samplesheet = samplesheet
            .groupTuple()  // Groups by sample_id: [sample_id, [R1_files...], [R2_files...]]
            .map { sample_id, r1_files, r2_files ->
                // Flatten R1 and R2 files into single list
                [sample_id, [r1_files, r2_files].flatten()]
            }
        
        // Pass all samples to CAT_FASTQ (it will decide merge vs symlink)
        CAT_FASTQ(samplesheet)
        samplesheet = CAT_FASTQ.out.reads
    }
    
    // --- PACBIO: Merge multiple flow cells/files ---
    if (data_type == 'pacbio') {
        samplesheet = samplesheet
            .groupTuple()  // Groups by sample_id: [sample_id, [reads_files...]]
            .map { sample_id, reads_files ->
                // Flatten in case of nested lists
                def flat_reads = reads_files.flatten()
                [sample_id, flat_reads]
            }
        
        // Pass to CAT_LONGREAD with explicit format
        def format = params.longread_format
        CAT_LONGREAD(samplesheet, format)
        samplesheet = CAT_LONGREAD.out.reads
    }
    
    // --- NANOPORE: Merge multiple flow cells/files ---
    if (data_type == 'nanopore') {
        samplesheet = samplesheet
            .groupTuple()  // Groups by sample_id: [sample_id, [reads_files...]]
            .map { sample_id, reads_files ->
                // Flatten in case of nested lists
                def flat_reads = reads_files.flatten()
                [sample_id, flat_reads]
            }
        
        // Pass to CAT_LONGREAD with explicit format
        def format = params.longread_format
        CAT_LONGREAD(samplesheet, format)
        samplesheet = CAT_LONGREAD.out.reads
    }
    
    // --- HYBRID: Merge Illumina lanes AND long read files ---
    if (data_type == 'hybrid') {
        samplesheet = samplesheet
            .groupTuple()  // Groups by sample_id: [sample_id, [R1s...], [R2s...], [longreads...]]
        
        // Branch into illumina and longread channels
        illumina_ch = samplesheet
            .map { sample_id, r1_files, r2_files, longread_files ->
                // Flatten R1 and R2 files for CAT_FASTQ
                [sample_id, [r1_files, r2_files].flatten()]
            }
        
        longread_ch = samplesheet
            .map { sample_id, r1_files, r2_files, longread_files ->
                // Flatten and deduplicate long read files
                def flat_longreads = longread_files.flatten().unique()
                [sample_id, flat_longreads]
            }
        
        // Merge illumina lanes
        CAT_FASTQ(illumina_ch)
        
        // Merge long read files with explicit format
        def format = params.longread_format
        CAT_LONGREAD(longread_ch, format)
        
        // Recombine: [sample_id, merged_R1, merged_R2, merged_longreads]
        samplesheet = CAT_FASTQ.out.reads
            .join(CAT_LONGREAD.out.reads)  // Join on sample_id
    }

    // =========================================
    // Route to appropriate workflow
    // =========================================
    
    // Collect software versions
    GET_SOFTWARE_VERSIONS()
    
    // Collect QC files for MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.versions.collect().ifEmpty([]))
    
    if (data_type == 'illumina') {
        IlluminaWorkflow(samplesheet)
        ch_multiqc_files = ch_multiqc_files.mix(IlluminaWorkflow.out.qc_files.collect().ifEmpty([]))
    } else if (data_type == 'pacbio') {
        PacBioWorkflow(samplesheet)
        ch_multiqc_files = ch_multiqc_files.mix(PacBioWorkflow.out.qc_files.collect().ifEmpty([]))
    } else if (data_type == 'nanopore') {
        NanoporeWorkflow(samplesheet)
        ch_multiqc_files = ch_multiqc_files.mix(NanoporeWorkflow.out.qc_files.collect().ifEmpty([]))
    } else if (data_type == 'hybrid') {
        HybridWorkflow(
            samplesheet,
            params.longread_type.toLowerCase(),
            params.longread_format.toLowerCase()
        )
        ch_multiqc_files = ch_multiqc_files.mix(HybridWorkflow.out.qc_files.collect().ifEmpty([]))
    }
    
    // Run MultiQC if we have QC files
    MULTIQC(ch_multiqc_files.collect())
}
