process coverStats {
    label "wfbacterialgenomes"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
    output:
        tuple val(meta), path("*fwd.regions.bed.gz"), emit: fwd
        tuple val(meta), path("*rev.regions.bed.gz"), emit: rev
        tuple val(meta), path("*total.regions.bed.gz"), emit: all

    """
    mosdepth -n --fast-mode --by 200 --flag 16 -t $task.cpus "${meta.alias}.fwd" align.bam
    mosdepth -n --fast-mode --by 200 --include-flag 16 -t $task.cpus "${meta.alias}.rev" align.bam
    mosdepth -n --fast-mode --by 200 -t $task.cpus "${meta.alias}.total" align.bam
    """
}
