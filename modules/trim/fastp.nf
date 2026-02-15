process FASTP {
    tag "fastp on $sample_id"
    label 'io_limited'
    publishDir "${params.outdir}/trimming/fastp/${sample_id}", mode: 'copy', overwrite: true

    publishDir "${params.outdir}/trimming/fastp/log_files/${sample_id}", 
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}.clean.1.fastq.gz"), path("${sample_id}.clean.2.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_fastp.html"), emit: html
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    path('.command.*')

    script:
    """
    fastp -w ${task.cpus} -q 25 \\
        --cut_front --cut_tail \\
        -l 20 --detect_adapter_for_pe \\
        -i ${fq1} -I ${fq2} \\
        -o ${sample_id}.clean.1.fastq.gz \\
        -O ${sample_id}.clean.2.fastq.gz \\
        -h ${sample_id}_fastp.html \\
        -j ${sample_id}_fastp.json \\
        ${params.fastp_opts ?: ''}
    """
}
