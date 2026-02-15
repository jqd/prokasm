process METASPADES {
    tag "metaSPAdes on $sample_id"
    label 'mem_veryhigh'
    publishDir "${params.outdir}/assembly/metaspades/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.scaffolds.fasta.gz'

    publishDir "${params.outdir}/assembly/metaspades/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(fq1), path(fq2)

    output:
    tuple val(sample_id), path("${sample_id}.scaffolds.fasta.gz"), emit: contigs
    path('.command.*')

    script:
    // Handle k-mer option
    def kmer_opt = params.metaspades_k == 'auto' ? '' : "-k ${params.metaspades_k}"
    """
    metaspades.py \\
      -1 $fq1 \\
      -2 $fq2 \\
      -o spades_out \\
      ${kmer_opt} \\
      -t ${task.cpus} \\
      ${params.spades_opts ?: ''}
    
    # Compress and rename output
    gzip -c spades_out/scaffolds.fasta > ${sample_id}.scaffolds.fasta.gz
    """
}
