process FLYE {
    tag "FLYE on $sample_id"
    label 'mem_veryhigh'
    publishDir "${params.outdir}/assembly/flye/${sample_id}", mode: 'copy', overwrite: true, pattern: "${sample_id}*"

    publishDir "${params.outdir}/assembly/flye/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.assembly.fasta.gz"), emit: contigs
    tuple val(sample_id), path("${sample_id}.assembly_graph.gfa.gz"), emit: gfa
    tuple val(sample_id), path("${sample_id}.assembly_graph.gv.gz"), emit: gv
    tuple val(sample_id), path("${sample_id}.assembly_info.txt"), emit: info
    tuple val(sample_id), path("${sample_id}.flye.log"), emit: log
    tuple val(sample_id), path("${sample_id}.params.json"), emit: params_json
    path('.command.*')

    script:
    def valid_modes = ["pacbio-raw", "pacbio-corr", "pacbio-hifi", "nano-raw", "nano-corr", "nano-hq"]
    if (!valid_modes.contains(params.flye_mode)) {
        error "Unrecognised mode to run Flye. Options: ${valid_modes.join(', ')}"
    }
    """
    flye --${params.flye_mode} ${reads} \\
        --out-dir flye_out \\
        --threads ${task.cpus} \\
        ${params.flye_opts ?: ''}

    # Compress and rename outputs
    gzip -c flye_out/assembly.fasta > ${sample_id}.assembly.fasta.gz
    gzip -c flye_out/assembly_graph.gfa > ${sample_id}.assembly_graph.gfa.gz
    gzip -c flye_out/assembly_graph.gv > ${sample_id}.assembly_graph.gv.gz
    mv flye_out/assembly_info.txt ${sample_id}.assembly_info.txt
    mv flye_out/flye.log ${sample_id}.flye.log
    mv flye_out/params.json ${sample_id}.params.json
    """
}
