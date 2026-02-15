process readStats {
    label "wfbacterialgenomes"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
    output:
        path "*readstats.txt", emit: stats
    """
    bamstats align.bam > "${meta.alias}.readstats.txt"
    if [[ \$(wc -l <"${meta.alias}.readstats.txt") -le 1 ]]; then
        echo "No alignments of reads to reference sequence found."
        exit 1
    fi
    """
}
