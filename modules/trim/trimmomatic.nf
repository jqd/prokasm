process TRIMMOMATIC {
    tag "trimmomatic on $sample_id"
    label 'io_limited'
    publishDir "${params.outdir}/trimming/trimmomatic/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.fastq.gz'
    publishDir "${params.outdir}/trimming/trimmomatic/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.{log,summary}'

    publishDir "${params.outdir}/trimming/trimmomatic/log_files/${sample_id}", 
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.paired_R1.fastq.gz"), path("${sample_id}.trimmed.paired_R2.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}.trim.log"), emit: log
    tuple val(sample_id), path("${sample_id}.summary"), emit: summary
    path('.command.*')

    script:
    def fq1_paired = "${sample_id}.trimmed.paired_R1.fastq.gz"
    def fq1_unpaired = "${sample_id}.trimmed.unpaired_R1.fastq.gz"
    def fq2_paired = "${sample_id}.trimmed.paired_R2.fastq.gz"
    def fq2_unpaired = "${sample_id}.trimmed.unpaired_R2.fastq.gz"
    
    // Default trimmomatic options if none provided
    def trim_opts = params.trimmomatic_opts ?: 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36'
    """
    trimmomatic PE -phred33 \\
        -threads ${task.cpus} \\
        -trimlog ${sample_id}.trim.log \\
        -summary ${sample_id}.summary \\
        ${fq1} \\
        ${fq2} \\
        ${fq1_paired} \\
        ${fq1_unpaired} \\
        ${fq2_paired} \\
        ${fq2_unpaired} \\
        ${trim_opts}
    """
}
