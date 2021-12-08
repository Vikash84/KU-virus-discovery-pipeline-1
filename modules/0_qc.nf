nextflow.enable.dsl=2

workflow qc_illumina {
    take:
        pe1
        pe2
    
    main:
        qc1 = fastqc1(pe1)
        qc2 = fastqc2(pe2)
        multiqc(qc1, qc2)
}

workflow qc_nanopore {
    take:
        fastq
    
    main:
        nanoplot(fastq)
}

process fastqc1 {
    conda '/home/molecularvirology/miniconda2/envs/vdp_srs'
    input:
        path fastq
    output:
        path "${fastq.simpleName}_fastqc.zip"
    """
    fastqc $fastq
    """
}

process fastqc2 {
    conda '/home/molecularvirology/miniconda2/envs/vdp_srs'
    input:
        path fastq
    output:
        path "${fastq.simpleName}_fastqc.zip"
    """
    fastqc $fastq 
    """
}

process multiqc {
    conda '/home/molecularvirology/miniconda2/envs/vdp_srs'
    publishDir "${params.outdir}/qc", mode: 'copy'
    input:
        path qc1
        path qc2
    output:
        path "*"
    """
    multiqc $qc1 $qc2
    """
}

process nanoplot {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    publishDir "${params.outdir}/qc", mode: 'copy'
    input:
        path fastq
    output:
        path "${params.prefix}_*"
    """
    NanoPlot --fasta $fastq -p ${params.prefix}_ -o ${params.prefix}
    mv ${params.prefix}/* .
    """
}
