process MEGAHIT {
    tag "megahit on $sample_id"
    label 'mem_medium'
    publishDir "${params.outdir}/assembly/megahit/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.contigs.fa.gz'

    publishDir "${params.outdir}/assembly/megahit/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.fa.gz"), emit: contigs
    path('.command.*')

    script:
    """
    megahit -t ${task.cpus} --force \\
        -1 ${fq1} \\
        -2 ${fq2} \\
        --out-prefix ${sample_id} \\
        -o megahit_out \\
        ${params.megahit_opts ?: ''}
    
    # Compress and rename output
    gzip -c megahit_out/${sample_id}.contigs.fa > ${sample_id}.contigs.fa.gz
    """
}
