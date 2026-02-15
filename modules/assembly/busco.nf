process BUSCO {
    tag "BUSCO on $sample_id"
    label 'io_net'
    publishDir "${params.outdir}/busco/", 
      mode: 'copy', 
      overwrite: true, 
      pattern: "${sample_id}_busco/*"

    publishDir "${params.outdir}/busco/log_files/${sample_id}",
      mode: 'copy',
      overwrite: true,
      pattern: '.command*',
      saveAs: { "${it}".replaceFirst(/^\./, "")}

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_busco/short_summary*"), emit: summary
    tuple val(sample_id), path("${sample_id}_busco/*.png"), optional: true, emit: plot
    path('.command.*')

    script:
    // Handle both gzipped and non-gzipped input
    def is_gzipped = contigs.toString().endsWith('.gz')
    
    // Build lineage option if specified
    def lineage_opt = params.busco_lineage ? "-l ${params.busco_lineage}" : "--auto-lineage-prok"
    """
    # Decompress if needed (BUSCO doesn't handle gzipped files well)
    if ${is_gzipped}; then
        gunzip -c ${contigs} > assembly.fasta
        INPUT_FILE="assembly.fasta"
    else
        INPUT_FILE="${contigs}"
    fi
    
    # Run BUSCO
    busco -i \${INPUT_FILE} \\
        -o ${sample_id}_busco \\
        -f \\
        -m ${params.busco_mode} \\
        ${lineage_opt} \\
        -c ${task.cpus} \\
        ${params.busco_opts ?: ''}
    
    # Generate plot (may fail if not all components are available)
    generate_plot.py -wd ${sample_id}_busco || true
    
    # Rename plot if it exists
    if [ -f "${sample_id}_busco/busco_figure.png" ]; then
        mv ${sample_id}_busco/busco_figure.png ${sample_id}_busco/${sample_id}_busco_figure.png
    fi
    """
}