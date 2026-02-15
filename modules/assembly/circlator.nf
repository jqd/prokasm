process CIRCLATOR {
    tag "CIRCLATOR on $sample_id"
    label 'io_net'
    publishDir "${params.outdir}/circlator/${sample_id}", mode: 'copy', overwrite: true, pattern: "${sample_id}*"

    publishDir "${params.outdir}/circlator/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_circlator_fixstart.fasta.gz"), emit: contigs
    tuple val(sample_id), path("${sample_id}_circlator_fixstart.log"), emit: log
    path('.command.*')

    script:
    // Handle both gzipped and non-gzipped input
    def is_gzipped = contigs.toString().endsWith('.gz')
    """
    # Decompress if needed
    if ${is_gzipped}; then
        gunzip -c ${contigs} > assembly.fasta
        INPUT_FILE="assembly.fasta"
    else
        INPUT_FILE="${contigs}"
    fi
    
    # Run circlator
    circlator fixstart \${INPUT_FILE} ${sample_id}_circlator_fixstart \\
        --min_id ${params.circlator_min_id} \\
        ${params.circlator_opts ?: ''}
    
    # Compress output for consistency with other tools
    gzip -c ${sample_id}_circlator_fixstart.fasta > ${sample_id}_circlator_fixstart.fasta.gz
    """
}
