process GET_SOFTWARE_VERSIONS {
    tag "Software versions"
    label 'low_cpu'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true

    output:
    path "software_versions.yml", emit: versions
    path "software_versions.txt", emit: versions_txt

    script:
    """
    echo "Pipeline: ProkASM v${workflow.manifest.version}" > software_versions.txt
    nextflow -v 2>&1 | grep "nextflow" >> software_versions.txt
    echo "" >> software_versions.txt
    echo "Software versions:" >> software_versions.txt
    
    # QC tools
    fastqc --version 2>&1 | grep "FastQC" >> software_versions.txt || echo "FastQC: Not available" >> software_versions.txt
    NanoPlot --version 2>&1 | grep "NanoPlot" >> software_versions.txt || echo "NanoPlot: Not available" >> software_versions.txt
    multiqc --version 2>&1 >> software_versions.txt || echo "MultiQC: Not available" >> software_versions.txt
    
    # Trimming tools
    fastp --version 2>&1 | head -n1 >> software_versions.txt || echo "fastp: Not available" >> software_versions.txt
    echo -n "Trimmomatic: " >> software_versions.txt && trimmomatic -version 2>&1 | tail -n 1 >> software_versions.txt || echo "Trimmomatic: Not available" >> software_versions.txt
    echo -n "Trim Galore: " >> software_versions.txt && trim_galore --version 2>&1 | grep "version" | sed 's/^[[:space:]]*version //' >> software_versions.txt || echo "Trim Galore: Not available" >> software_versions.txt
    
    # Assemblers
    megahit --version 2>&1 | head -n1 >> software_versions.txt || echo "MEGAHIT: Not available" >> software_versions.txt
    spades.py --version 2>&1 | head -n1 >> software_versions.txt || echo "SPAdes: Not available" >> software_versions.txt
    echo -n "Hifiasm: " >> software_versions.txt && hifiasm --version 2>&1 | head -n1 >> software_versions.txt || echo "Hifiasm: Not available" >> software_versions.txt
    echo -n "Flye: " >> software_versions.txt && flye --version 2>&1 >> software_versions.txt || echo "Flye: Not available" >> software_versions.txt
    unicycler --version 2>&1 >> software_versions.txt || echo "Unicycler: Not available" >> software_versions.txt
    
    # Post-assembly tools
    dnaapler --version 2>&1 >> software_versions.txt || echo "dnaapler: Not available" >> software_versions.txt
    echo -n "Circlator: " >> software_versions.txt && circlator version 2>&1 >> software_versions.txt || echo "Circlator: Not available" >> software_versions.txt
    busco --version 2>&1 >> software_versions.txt || echo "BUSCO: Not available" >> software_versions.txt
    
    # Alignment/mapping tools
    echo -n "minimap2: " >> software_versions.txt && minimap2 --version 2>&1 >> software_versions.txt || echo "minimap2: Not available" >> software_versions.txt
    samtools --version 2>&1 | head -n1 >> software_versions.txt || echo "samtools: Not available" >> software_versions.txt
    
    # Utility tools
    bamtools --version 2>&1 | grep "bamtools" >> software_versions.txt || echo "bamtools: Not available" >> software_versions.txt
    
    # Create YAML version for MultiQC
    cat <<EOF > software_versions.yml
    pipeline:
      prokasm: "${workflow.manifest.version}"
      nextflow_info: "\$(nextflow info 2>&1 | tr '\n' ' ')"
    software:
EOF
    
    # Add each tool to YAML (simplified)
    grep -E "^[A-Za-z]" software_versions.txt | grep -v "Pipeline:" | grep -v "Nextflow:" | grep -v "Software versions:" | awk '{print "  " tolower(\$1) ": " \$0}' >> software_versions.yml || true
    """
}
