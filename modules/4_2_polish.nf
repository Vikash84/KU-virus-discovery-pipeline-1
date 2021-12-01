nextflow.enable.dsl=2

workflow polish {
    take:
        contigs
        fastq
    emit:
        polished_contigs
    main:
        corrected_contigs = racon(contigs, fastq)
        polished_contigs = medaka(corrected_contigs, fastq)
}

process racon {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
        path contigs
        path fastq
    output:
        path "${params.prefix}.polished_contigs.fasta"
    """
    echo "Contigs correction with Racon"
    echo "Mapping reads onto contigs using Minimap2"
    minimap2 -x ava-ont $contigs $fastq \
        -t 24 > ${params.prefix}.paf
    echo "Correcting contigs using Racon"
    racon -t 24 $fastq ${params.prefix}.paf $contigs > ${params.prefix}.polished_contigs.fasta
    """
}

process medaka {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
        path contigs
        path fastq
    output:
        path "${params.prefix}/${params.prefix}_consensus.fasta"
    """
    echo "Contigs polishing with Medaka"
    medaka_consensus -i $fastq -d $contigs \
         -o ${params.prefix} \
         -t 24 -m r941_min_high_g303
    mv ${params.prefix}/consensus.fasta ${params.prefix}/${params.prefix}_consensus.fasta
    """
}