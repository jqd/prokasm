process MERGE_FASTQ {
    tag "MERGE_FASTQ on $sample_id"
    label = 'io_net'
    publishDir "${params.outdir}/merged_fastq/${sample_id}", mode: 'copy', overwrite: true

    publishDir "${params.outdir}/merged_fastq/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    tuple val(sample_id), path("${sample_id}_merged.fastq.gz"), emit: merged_fastq
    path('.command.*')

    script:
    """
    cat ${fastq_files.join(' ')} > ${sample_id}_merged.fastq.gz
    """
}
