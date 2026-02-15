process CAT_FASTQ {
    tag "concatenate on $sample_id"
    label 'io_limited'
    publishDir "${params.outdir}/concatenated/${sample_id}", mode: 'copy'
    publishDir "${params.outdir}/concatenated/log_files/${sample_id}", pattern: '.command*', mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1.merged.fastq.gz"), path("${sample_id}_2.merged.fastq.gz"), emit: reads
    
    script:
    def readList = reads instanceof List ? reads.collect { it.toString() } : [reads.toString()]
    
    // Detect if files are gzipped
    def isGzipped = readList[0].endsWith('.gz')
    
    if (readList.size() == 2) {
        // Single lane - create symlinks (ensure gzipped output)
        if (isGzipped) {
            """
            ln -s ${reads[0]} ${sample_id}_1.merged.fastq.gz
            ln -s ${reads[1]} ${sample_id}_2.merged.fastq.gz
            """
        } else {
            """
            gzip -c ${reads[0]} > ${sample_id}_1.merged.fastq.gz
            gzip -c ${reads[1]} > ${sample_id}_2.merged.fastq.gz
            """
        }
    } else if (readList.size() > 2) {
        // Multiple lanes - concatenate files
        // Separate R1 and R2 files (alternate indexing: 0,2,4... = R1; 1,3,5... = R2)
        def read1 = []
        def read2 = []
        readList.eachWithIndex { v, ix -> (ix & 1 ? read2 : read1) << v }
        
        if (isGzipped) {
            // For gzipped files, cat works correctly (concatenates gzip streams)
            """
            cat ${read1.join(' ')} > ${sample_id}_1.merged.fastq.gz
            cat ${read2.join(' ')} > ${sample_id}_2.merged.fastq.gz
            """
        } else {
            // For uncompressed files, concatenate and compress
            """
            cat ${read1.join(' ')} | gzip -c > ${sample_id}_1.merged.fastq.gz
            cat ${read2.join(' ')} | gzip -c > ${sample_id}_2.merged.fastq.gz
            """
        }
    } else {
        error("Could not find any FASTQ file pairs to concatenate in the process input")
    }
    
}
