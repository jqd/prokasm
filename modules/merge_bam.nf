process MERGE_BAM {
    tag "merge BAM on $sample_id"
    label 'io_net'
    publishDir "${params.outdir}/merged_bam/${sample_id}", mode: 'copy'
    publishDir "${params.outdir}/merged_bam/log_files/${sample_id}", pattern: '.command*', mode: 'copy'

    input:
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path("${sample_id}.merged.bam"), emit: bam
    
    script:
    def bamList = bams instanceof List ? bams.collect { it.toString() }.join(' ') : bams.toString()
    if (bams instanceof List && bams.size() > 1) {
        """
        samtools merge -@ ${task.cpus} ${sample_id}.merged.bam ${bamList}
        """
    } else {
        """
        ln -s ${bamList} ${sample_id}.merged.bam
        """
    }
}
