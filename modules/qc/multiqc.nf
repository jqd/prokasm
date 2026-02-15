process MULTIQC {
    tag "MultiQC aggregation"
    label 'io_limited'
    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data
    path('.command.*')

    script:
    """
    multiqc . \
        --force \
        --title "ProkASM Assembly Report" \
        --filename multiqc_report.html \
        --config /dev/null
    """
}
