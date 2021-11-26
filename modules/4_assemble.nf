nextflow.enable.dsl=2

workflow assemble_illumina {
    take:
        fastq
    emit:
        contigs
    main:
        contigs = spades(fastq)
        filtered_contigs = filterContigs_illumina(contigs)
        contigSummary(filtered_contigs)
        
}

workflow assemble_nanopore {
    take:
        fastq
    emit:
        contigs
    main:
        contigs = megahit(fastq)
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
    megahit -r $fastq -o $params.prefix --out-prefix $params.prefix -t 12 
    mv ${params.prefix}/${params.prefix}.contigs.fa ${params.prefix}.megahit.contigs.fa
    """
}

process canu {
    errorStrategy { 'ignore' }
    publishDir "${params.outdir}/assembly", mode: 'copy'
    
    input:
        path fastq
    output:
        path("${params.prefix}.canu.contigs.fa") optional true
    """
    echo "De novo assembly with Canu"
    canu -p $params.prefix -d $params.prefix -nanopore $fastq \
    genomeSize=5m minReadLength=300 minOverlapLength=50 maxThreads=12 minInpuCoverage=0 StopOnCoverage=0 \
    corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 \
    redMemory=32 oeaMemory=32 batMemory=200
    """
}