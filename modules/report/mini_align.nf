process MINI_ALIGN {
    tag "MINI_ALIGN on $sample_id"
    label 'multithread'
    publishDir "${params.outdir}/alignment/${sample_id}", mode: 'copy', overwrite: true

    publishDir "${params.outdir}/alignment/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(reads), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.stats.txt"), emit: stats
    path('.command.*')

    script:
    // Handle gzipped contigs
    def is_gzipped = contigs.toString().endsWith('.gz')
    """
    # Decompress contigs if needed
    if ${is_gzipped}; then
        gunzip -c ${contigs} > reference.fasta
    else
        cp ${contigs} reference.fasta
    fi
    
    # Index reference
    samtools faidx reference.fasta
    
    # Align reads to contigs using minimap2 and sort
    minimap2 -ax map-ont -t ${task.cpus} reference.fasta ${reads} | \\
        samtools view -bS - | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    
    # Index BAM
    samtools index ${sample_id}.sorted.bam
    
    # Generate alignment stats
    samtools stats ${sample_id}.sorted.bam > ${sample_id}.stats.txt
    """
}
