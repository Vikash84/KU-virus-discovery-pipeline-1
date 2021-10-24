nextflow.enable.dsl=2


workflow reference_map_illumina {
    take:
        fastq
    main:
        references = Channel.fromPath(params.ref_vir_list)
                            .splitText()
                            .map{file(it)}

        fastq_pair = split_fastq_tuple(fastq)
        bams = each_ref_mapping_illumina(fastq_pair[0].first(), fastq_pair[1].first(), references).filter{it.size()>1000}
        map_info = mapping_summary(bams)
        collection = map_info[0].collectFile(name: "${params.prefix}_temp_colleciton.txt")
        summary_collection = add_header(collection)
        filter(summary_collection)
}

workflow reference_map_nanopore {
    take:
        filtered
        ref_list

    main:
        Channel.fromPath(ref_list)
            .splitText().map { it -> it.trim() }
            .map { file(it) }
            .set { references }

        fastq_ref = filtered.combine(references)
        bams = each_ref_mapping_nanopore(fastq_ref).filter{it.size()>1000}
        map_info = mapping_summary(bams)
        collection = map_info[0].collectFile(name: "${params.prefix}_temp_colleciton.txt")
        summary_collection = add_header(collection)
        filter(summary_collection)
}

process generate_ref_vir_channel {
    input:
        tuple path(pe1), path(pe2)
    output:
        Channel.fromPath('/tmp/Viral_sequences_DB/*/*fas')
    """
    """
}

process split_fastq_tuple {
    input:
        tuple path(pe1), path(pe2)
    output:
        path("1.fq")
        path("2.fq")
    """
        cat $pe1 > 1.fq
        cat $pe2 > 2.fq
    """
}

process each_ref_mapping_illumina {
    errorStrategy 'ignore'
    stageInMode "link"

    input:
        path(pe1)
        path(pe2)
        path(ref_fa)
    output:
        path "${params.prefix}_${ref_fa.simpleName}.bam" optional true

    """
    bowtie2-build $ref_fa ${ref_fa.baseName}_index
    bowtie2 -1 $pe1 -2 $pe2 -x ${ref_fa.baseName}_index | samtools view -Sb - | samtools sort - | samtools view -b -F 4 -o "${params.prefix}_${ref_fa.simpleName}.bam"
    """
}

process each_ref_mapping_nanopore {

    input:
        tuple path(fastq), path(ref_fa)
    output:
        path "${params.prefix}_${ref_fa.simpleName}.bam" optional true

    """
    minimap2 -ax map-ont $ref_fa $fastq | samtools view -Sb - | samtools sort - | samtools view -b -F 4 -o "${params.prefix}_${ref_fa.simpleName}.bam"
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
    ~/bamcov -H $bam > each_mapping/"${bam.simpleName}.txt"
    qualimap bamqc -bam $bam -outdir qualimap_results
    mv qualimap_results/images_qualimapReport/genome_coverage_across_reference.png "${bam.simpleName}.png"
    """
}

process filter {
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
        path collection
    output:
        path "${params.prefix}.filtered_reference_mapping_collection.txt"
    
    """
    python ~/scripts/2_reference_map_result_filter.py --input $collection --output ${params.prefix}.filtered_reference_mapping_collection.txt --min_avg_cov 1.0
    """
}

process add_header {
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
        path collection
    output:
        path "${params.prefix}.reference_mapping_summary.txt"
    """
    cat ~/headers/bamcov_header $collection > ${params.prefix}.reference_mapping_summary.txt
    """
}
