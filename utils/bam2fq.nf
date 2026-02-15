process BAM2FQ {
    tag "$sample_id"
    label 'io_limited'
    publishDir "${params.outdir}/bam2fq/${sample_id}", mode: 'copy', overwrite: true

    publishDir "${params.outdir}/bam2fq/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.fastq.gz"), emit: fastq
    path('.command.*')

    script:
    """
    # Index BAM file first
    pbindex ${bam}
    
    # Convert to FASTQ (bam2fastq from pbtk)
    bam2fastq -c 1 -o ${sample_id} ${bam}
    """
}
