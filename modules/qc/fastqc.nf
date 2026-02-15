process FASTQC {
    tag "$sample_id"
    label 'io_limited'
    publishDir "${params.outdir}/qc/fastqc/${sample_id}", mode: 'copy', overwrite: true

    publishDir "${params.outdir}/qc/fastqc/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip
    path('.command.*')

    script:
    """
    fastqc -t ${task.cpus} ${read1} ${read2} -o . ${params.fastqc_opts ?: ''}
    """
}
