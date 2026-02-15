process TRIMGALORE {
    tag "trim_galore on $sample_id"
    label 'io_limited'
    publishDir "${params.outdir}/trimming/trimgalore/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.fq.gz'
    publishDir "${params.outdir}/trimming/trimgalore/${sample_id}/fastqc", mode: 'copy', overwrite: true, pattern: '*.{txt,zip,html}'

    publishDir "${params.outdir}/trimming/trimgalore/log_files/${sample_id}", 
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: reads
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip
    tuple val(sample_id), path("*_trimming_report.txt"), emit: report
    path('.command.*')

    script:
    """
    trim_galore \\
        --cores ${task.cpus} \\
        --quality ${params.trimgalore_quality} \\
        --stringency ${params.trimgalore_stringency} \\
        --paired \\
        --gzip \\
        --fastqc \\
        ${params.trimgalore_opts ?: ''} \\
        ${fq1} \\
        ${fq2}
    """
}
