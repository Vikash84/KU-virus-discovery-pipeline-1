nextflow.enable.dsl=2


workflow reference_mapping_illumina {
    take:
        fastq_pair
    main:
        bams = mapping_illumina(fastq_pair).flatten().filter{it.size()>1000}
        mapping_summaries = mapping_summary(bams)
        collection = mapping_summaries[0].collectFile(name: "${params.prefix}_temp_colleciton.txt")
        add_header_filter(collection)
}

workflow reference_mapping_nanopore {
    take:
        fastq
    main:
        bams = mapping_nanopore(fastq).flatten().filter{it.size()>1000}
        mapping_summaries = mapping_summary(bams)
        collection = mapping_summaries[0].collectFile(name: "${params.prefix}_temp_colleciton.txt")
        add_header_filter(collection)
}

process mapping_illumina {
    stageInMode "link"

    input:
        tuple path(pe1), path(pe2)
    output:
        file '*bam' optional true
    """
    cat ${params.reference_list_path} | while read ref
    do
    base=\$(basename \${ref})
    simple=\${base%%.*}
    minimap2 -ax sr \${ref} $pe1 $pe2 | samtools view -Sb | samtools sort - | samtools view -b -F 4 -o "${params.prefix}_\${simple}.bam"
    done
    """
}

process mapping_nanopore {
    stageInMode "link"

    input:
        path(fastq)
    output:
        file "*.bam" optional true

    """
    cat ${params.reference_list_path} | while read ref
    do
    base=\$(basename \${ref})
    simple=\${base%%.*}
    minimap2 -ax map-ont \${ref} $fastq | samtools view -Sb - | samtools sort - | samtools view -b -F 4 -o "${params.prefix}_\${simple}.bam"
    done
    """
}

process mapping_summary {
    publishDir "${params.outdir}/mapping/details", mode: 'copy'
    errorStrategy 'ignore'
    stageInMode 'link'

    input:
        path bam
    output:
        path "${bam.simpleName}.txt"
        path "${bam.simpleName}/*"
    """
    ${params.bamcov_path}/bamcov -H $bam > "${bam.simpleName}.txt"
    qualimap bamqc -bam $bam -outdir ${bam.simpleName}
    """
}

process add_header_filter {
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
        path collection
    output:
        path "${params.prefix}.reference_mapping.txt"
    
    """
    cat ${params.nextflow_script_path}/headers/bamcov_header $collection > ${params.prefix}.tmp
    python ${params.nextflow_script_path}/scripts/2_reference_mapping_result_filter.py --input ${params.prefix}.tmp --output ${params.prefix}.reference_mapping.txt --min_avg_cov ${params.reference_mapping_minimum_avg_coverage}
    """
}
