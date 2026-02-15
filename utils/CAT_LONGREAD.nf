/**
 * CAT_LONGREAD - Concatenate or symlink long read files
 * 
 * Handles merging of multiple long read files per sample
 * Uses provided format parameter to determine handling (BAM vs FASTQ)
 */
process CAT_LONGREAD {
    tag "concatenate longreads on $sample_id"
    label 'io_limited'
    publishDir "${params.outdir}/concatenated_longreads/${sample_id}", mode: 'copy'
    publishDir "${params.outdir}/concatenated_longreads/log_files/${sample_id}", pattern: '.command*', mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val format

    output:
    tuple val(sample_id), path("${sample_id}.merged.*"), emit: reads
    
    script:
    def readList = reads instanceof List ? reads.collect { it.toString() } : [reads.toString()]
    
    // Use explicit format parameter
    def isBam = format == 'bam'
    
    // Check compression for FASTQ handling
    def firstFile = readList[0].toLowerCase()
    def isGzipped = firstFile.endsWith('.gz')
    
    if (readList.size() == 1) {
        // Single file - create symlink (or compress if needed)
        if (isBam) {
            """
            ln -s ${reads} ${sample_id}.merged.bam
            """
        } else {
            // FASTQ handling - standardize to .fastq.gz
            if (isGzipped) {
                """
                ln -s ${reads} ${sample_id}.merged.fastq.gz
                """
            } else {
                """
                gzip -c ${reads} > ${sample_id}.merged.fastq.gz
                """
            }
        }
    } else {
        // Multiple files - concatenate
        if (isBam) {
            // For BAM files, use samtools merge
            """
            samtools merge -@ ${task.cpus} ${sample_id}.merged.bam ${reads}
            """
        } else {
            // For FASTQ files
            if (isGzipped) {
                // For gzipped FASTQ, cat works (multi-stream gzip)
                """
                cat ${readList.join(' ')} > ${sample_id}.merged.fastq.gz
                """
            } else {
                // For uncompressed FASTQ, concatenate and compress
                """
                cat ${readList.join(' ')} | gzip -c > ${sample_id}.merged.fastq.gz
                """
            }
        }
    }
}
