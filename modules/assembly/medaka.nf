process FLYE {
    tag "FLYE on $sample_id"
    label = 'io_net'
    publishDir "${params.outdir}/flye/${sample_id}", mode: 'copy', overwrite: true, pattern: "${sample_id}*"

    publishDir "${params.outdir}/flye/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}
    input:
    tuple val(sample_id), path(read1)

    output:
    tuple val(sample_id), path("*.fasta.gz"), emit: contigs
    tuple val(sample_id), path("*.gfa.gz")  , emit: gfa
    tuple val(sample_id), path("*.gv.gz")   , emit: gv
    tuple val(sample_id), path("*.txt")     , emit: txt
    tuple val(sample_id), path("*.log")     , emit: log
    tuple val(sample_id), path("*.json")    , emit: json
    path('.command.*')

    script:
    """
    flye ${params.flye_mode} $read1 --out-dir .  --threads ${task.cpus} ${params.flye_opts}

    gzip -c assembly.fasta > ${sample_id}.assembly.fasta.gz
    gzip -c assembly_graph.gfa > ${sample_id}.assembly_graph.gfa.gz
    gzip -c assembly_graph.gv > ${sample_id}.assembly_graph.gv.gz
    mv assembly_info.txt ${sample_id}.assembly_info.txt
    mv flye.log ${sample_id}.flye.log
    mv params.json ${sample_id}.params.json

    """
}
