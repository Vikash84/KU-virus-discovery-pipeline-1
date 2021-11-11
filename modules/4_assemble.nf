nextflow.enable.dsl=2

workflow assemble_illumina {
    take:
        fastq
    emit:
        contigs
    main:
        spades_contigs = spades(fastq)
        contigs = spades_contigs
        filtered_contigs = filterContigs_illumina(contigs)
        contigSummary(filtered_contigs)
        
}

workflow assemble_nanopore {
    take:
        fastq
    emit:
        contigs
    main:
        megahit_contigs = megahit(fastq)
        contigs = megahit_contigs
        filtered_contigs = filterContigs_nanopore(contigs)
        contigSummary(filtered_contigs)
}

process filterContigs_illumina {
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
        path contigs
    output:
        path("${params.prefix}.filtered_contigs.fasta")
    """
    reformat.sh in=$contigs out="${params.prefix}.filtered_contigs.fasta" minlength=${params.illumina_min_contig_length}
    """
}

process filterContigs_nanopore {
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
        path contigs
    output:
        path("${params.prefix}.filtered_contigs.fasta")
    """
    reformat.sh in=$contigs out="${params.prefix}.filtered_contigs.fasta" minlength=${params.nanopore_min_contig_length}
    """
}

process contigSummary {
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
        path contigs
    output:
        path "${params.prefix}.contig_summary.txt"
        path "${params.prefix}.contigs_length_histogram.png"
    """
    python ${params.nextflow_script_path}/scripts/4_contig_summary_statistics.py $params.prefix $contigs
    """
}

process spades{
    errorStrategy { 'ignore' }
    input:
        tuple path(pe1), path(pe2)
    output:
        path("${params.prefix}.spades.contigs.fasta") optional true
    """
    echo "De novo assembly with Spades"
    echo "input file 1: $pe1"
    echo "input file 2: $pe2"
    spades.py --pe1-1 $pe1 --pe1-2 $pe2 -o ${params.prefix} --meta -t 12
    mv "${params.prefix}/contigs.fasta" "${params.prefix}.spades.contigs.fasta"
    """
}

process megahit {
    errorStrategy { 'ignore' }
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
        path fastq
    output:
        path("${params.prefix}.megahit.contigs.fa") optional true
    """
    echo "De novo assembly with Megahit"
    echo "input file: $fastq"
    megahit -r $fastq -o $params.prefix --out-prefix $params.prefix -t 12 
    mv ${params.prefix}/${params.prefix}.contigs.fa ${params.prefix}.megahit.contigs.fa
    """
}
