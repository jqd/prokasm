nextflow -v 2>&1 | grep "nextflow" >> software_versions.txt
    echo "" >> software_versions.txt
    echo "Software versions:" >> software_versions.txt
    
    # QC tools
    fastqc --version 2>&1 | grep "FastQC" >> software_versions.txt || echo "FastQC: Not available" >> software_versions.txt
    NanoPlot --version 2>&1 | grep "NanoPlot" >> software_versions.txt || echo "NanoPlot: Not available" >> software_versions.txt
    multiqc --version 2>&1 >> software_versions.txt || echo "MultiQC: Not available" >> software_versions.txt
    
    # Trimming tools
    fastp --version 2>&1 | head -n1 >> software_versions.txt || echo "fastp: Not available" >> software_versions.txt
    trimmomatic -version 2>&1 | head -n1 >> software_versions.txt || echo "Trimmomatic: Not available" >> software_versions.txt
    trim_galore --version 2>&1 | head -n1 >> software_versions.txt || echo "Trim Galore: Not available" >> software_versions.txt
    
    # Assemblers
    megahit --version 2>&1 | head -n1 >> software_versions.txt || echo "MEGAHIT: Not available" >> software_versions.txt
    spades.py --version 2>&1 | head -n1 >> software_versions.txt || echo "SPAdes: Not available" >> software_versions.txt
    hifiasm --version 2>&1 | head -n1 >> software_versions.txt || echo "Hifiasm: Not available" >> software_versions.txt
    flye --version 2>&1 >> software_versions.txt || echo "Flye: Not available" >> software_versions.txt
    unicycler --version 2>&1 >> software_versions.txt || echo "Unicycler: Not available" >> software_versions.txt
    
    # Post-assembly tools
    dnaapler --version 2>&1 >> software_versions.txt || echo "dnaapler: Not available" >> software_versions.txt
    circlator version 2>&1 >> software_versions.txt || echo "Circlator: Not available" >> software_versions.txt
    busco --version 2>&1 >> software_versions.txt || echo "BUSCO: Not available" >> software_versions.txt
    
    # Alignment/mapping tools
    minimap2 --version 2>&1 >> software_versions.txt || echo "minimap2: Not available" >> software_versions.txt
    samtools --version 2>&1 | head -n1 >> software_versions.txt || echo "samtools: Not available" >> software_versions.txt
    
    # Utility tools
    bamtools --version 2>&1 | head -n1 >> software_versions.txt || echo "bamtools: Not available" >> software_versions.txt
