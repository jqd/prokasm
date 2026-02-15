process DNAAPLER {
    tag "DNAAPLER on $sample_id"
    label 'io_net'
    publishDir "${params.outdir}/dnaapler/${sample_id}", mode: 'copy', overwrite: true, pattern: "${sample_id}*"

    publishDir "${params.outdir}/dnaapler/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_dnaapler_reoriented.fasta.gz"), emit: contigs
    path('.command.*')

    script:
    // Handle both gzipped and non-gzipped input
    def is_gzipped = contigs.toString().endsWith('.gz')
    """
    # Decompress if needed
    if ${is_gzipped}; then
        gunzip -c ${contigs} > assembly.fasta
    else
        cp ${contigs} assembly.fasta
    fi
    
    # Run dnaapler
    dnaapler all \\
        --input assembly.fasta \\
        --output dnaapler_output \\
        --prefix ${sample_id}_dnaapler \\
        --threads ${task.cpus} \\
        ${params.dnaapler_opts ?: ''}
    
    # Compress output
    gzip -c "dnaapler_output/${sample_id}_dnaapler_reoriented.fasta" > "${sample_id}_dnaapler_reoriented.fasta.gz"
    """
}
