process HIFIASM {
    tag "HIFIASM on $sample_id"
    label 'mem_veryhigh'
    publishDir "${params.outdir}/assembly/hifiasm/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.fa.gz'
    publishDir "${params.outdir}/assembly/hifiasm/${sample_id}", mode: 'copy', overwrite: true, pattern: '*.gfa.gz'

    publishDir "${params.outdir}/assembly/hifiasm/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.asm.bp.p_ctg.fa.gz"), emit: contigs
    tuple val(sample_id), path("${sample_id}.asm.bp.p_ctg.gfa.gz"), emit: gfa
    path('.command.*')

    script:
    """
    hifiasm -o ${sample_id}.asm \\
      -t ${task.cpus} \\
      -k ${params.hifiasm_k} \\
      -a ${params.hifiasm_a} \\
      ${params.hifiasm_opts ?: ''} \\
      ${reads}
    
    # Convert GFA to FASTA
    awk '/^S/{print ">"\$2;print \$3}' ${sample_id}.asm.bp.p_ctg.gfa > ${sample_id}.asm.bp.p_ctg.fa
    
    # Compress outputs
    gzip -c ${sample_id}.asm.bp.p_ctg.fa > ${sample_id}.asm.bp.p_ctg.fa.gz
    gzip -c ${sample_id}.asm.bp.p_ctg.gfa > ${sample_id}.asm.bp.p_ctg.gfa.gz
    """
}
