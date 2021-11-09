nextflow.enable.dsl=2


workflow reference_map_illumina {
    take:
        fastq
    main:
        bams = ref_mapping_illumina(fastq).filter{it.size()>1000}
        map_info = mapping_summary(bams)
        collection = map_info[0].collectFile(name: "${params.prefix}_temp_colleciton.txt")
        summary_collection = add_header(collection)
        add_header_filter(summary_collection)
}

workflow reference_map_nanopore {
    take:
        fastq

    main:
        bams = ref_mapping_nanopore(fastq).filter{it.size()>1000}
        map_info = mapping_summary(bams)
        collection = map_info[0].collectFile(name: "${params.prefix}_temp_colleciton.txt")
        summary_collection = add_header(collection)
        add_header_filter(summary_collection)
}

process ref_mapping_illumina {
    errorStrategy 'ignore'
    stageInMode "link"

    input:
        tuple path(pe1), path(pe2)
    output:
        file '*bam'
    """
    cat ${params.reference_list_path} | while read ref
    do
    base=\$(basename \${ref})
    simple=\${base%%.*}
    bowtie2 -1 $pe1 -2 $pe2 -x \${ref} | samtools view -Sb | samtools sort - | samtools view -b -F 4 -o "${params.prefix}_\${simple}.bam"
    done
    """
}

process ref_mapping_nanopore {

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
    errorStrategy 'ignore'
    stageInMode 'link'

    input:
        path bam
    output:
        path "each_mapping/${bam.simpleName}.txt"
        path "${bam.simpleName}.png" optional true
    """
    mkdir each_mapping
    bamcov -H $bam > each_mapping/"${bam.simpleName}.txt"
    qualimap bamqc -bam $bam -outdir qualimap_results
    mv qualimap_results/images_qualimapReport/genome_coverage_across_reference.png "${bam.simpleName}.png"
    """
}

process add_header_filter {
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
        path collection
    output:
        path "${params.prefix}.reference_mapping_collection.txt"
    
    """
    cat ~/headers/bamcov_header $collection > ${params.prefix}.tmp
    python ~/scripts/2_reference_map_result_filter.py --input ${params.prefix}.tmp --output ${params.prefix}.reference_mapping_collection.txt --min_avg_cov ${params.reference_mapping_minimum_avg_coverage}
    """
}
