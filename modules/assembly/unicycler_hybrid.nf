process UNICYCLER_HYBRID {
    tag "unicycler hybrid on $sample_id"
    label 'mem_veryhigh'
    publishDir "${params.outdir}/assembly/unicycler_hybrid/${sample_id}", mode: 'copy', overwrite: true

    publishDir "${params.outdir}/assembly/unicycler_hybrid/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(read1), path(read2), path(longreads)

    output:
    tuple val(sample_id), path("${sample_id}.assembly.fasta.gz"), emit: contigs
    tuple val(sample_id), path("${sample_id}.assembly.gfa.gz"), emit: gfa
    path("${sample_id}.unicycler.log")
    path('.command.*')

    script:
    """
    unicycler -t ${task.cpus} \\
        --mode ${params.unicycler_mode ?: 'normal'} \\
        --min_fasta_length ${params.unicycler_min_fasta_length ?: 100} \\
        -1 ${read1} \\
        -2 ${read2} \\
        -l ${longreads} \\
        -o unicycler_out \\
        ${params.unicycler_opts ?: ''}
    
    # Compress and rename outputs
    gzip -c unicycler_out/assembly.fasta > ${sample_id}.assembly.fasta.gz
    gzip -c unicycler_out/assembly.gfa > ${sample_id}.assembly.gfa.gz
    cp unicycler_out/unicycler.log ${sample_id}.unicycler.log
    """
}
