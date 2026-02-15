process NANOPLOT {
    tag "$sample_id"
    label 'io_limited'
    publishDir "${params.outdir}/qc/nanoplot/${sample_id}", mode: 'copy', overwrite: true

    publishDir "${params.outdir}/qc/nanoplot/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_NanoPlot-report.html"), emit: report
    tuple val(sample_id), path("${sample_id}_NanoStats.txt"), emit: stats
    tuple val(sample_id), path("*.png"), emit: plots
    path('.command.*')

    script:
    """
    NanoPlot --fastq ${reads} \\
        --threads ${task.cpus} \\
        --prefix ${sample_id}_ \\
        -o . \\
        ${params.nanoplot_opts ?: ''}
    """
}
