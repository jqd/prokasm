process UNICYCLER {
    tag "unicycler on $sample_id"
    label = 'mem_medium'
    publishDir "${params.outdir}/unicycler/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.contigs.fa'

    publishDir "${params.outdir}/unicycler/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: ['.command*', '*.log']
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}.assembly.fasta"), emit: contigs
    path('.command.*')
    path('unicycler.log')

    script:
    """
    unicycler -t ${task.cpus} \
        --mode ${params.unicycler_mode ?: 'normal'} \
        --min_fasta_length ${params.unicycler_min_fasta_length ?: 100} \
        -1 ${fq1} \
        -2 ${fq2} \
        -o . \
        ${params.unicycler_opts ?: ''}
    cp assembly.fasta  ${sample_id}.assembly.fasta
    """
}
