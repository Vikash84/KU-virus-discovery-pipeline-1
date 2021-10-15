nextflow.enable.dsl=2

workflow filter_illumina {
    take:
        pe1
        pe2
    emit:
        filtered
    main:
        decompressed = decompress_fastq(pe1, pe2)
        dedup = deduplication(decompressed) 
        filtered = trimmomatic(dedup)
        
        if ( params.host )
            filtered = hostfilter_illumina(filtered)
}

workflow filter_nanopore {
    take:
        fastq
    emit:
        filtered
    main:
        filtered = nanofilt(fastq)
        
        if ( params.host )
            filtered = hostfilter_nanopore(filtered)
}

process decompress_fastq {
    stageInMode "copy"

    input:
        path pe1
        path pe2
    output:
        tuple path("${pe1.baseName}"), path("${pe2.baseName}")
    
    script:
    """
    if [[ $pe1 == *.gz ]]; then
        gunzip $pe1
    fi

    if [[ $pe2 == *.gz ]]; then
        gunzip $pe2
    fi
    """
    
}

process deduplication {
    label "containerPrinSeq"
    conda "/home/molecularvirology/miniconda2/envs/vdp_srs"

    input:
        tuple path(pe1), path(pe2)
    output:
        tuple path("${params.prefix}_dedup_1.fastq"), path("${params.prefix}_dedup_2.fastq")
    script:
    """
    prinseq-lite.pl -fastq $pe1 -fastq2 $pe2 -derep 12 -out_good ${params.prefix}_dedup -out_bad null
    """
}

process trimmomatic {
    label "containerTrimmomatic"
    conda "/home/molecularvirology/miniconda2/envs/vdp_srs"
    input:
        tuple path(pe1), path(pe2)
    output:
        tuple path("${params.prefix}_trimmed_1.fq"), path ("${params.prefix}_trimmed_2.fq")
    script:
        
    """
    echo "Filter low quality reads"
    echo "sample name: ${params.prefix}"
    echo "input file 1: $pe1"
    echo "input file 2: $pe2"
    
    trimmomatic PE -phred33 $pe1 $pe2 ${params.prefix}_trimmed_1.fq unpaired.fq ${params.prefix}_trimmed_2.fq unpaired.fq ILLUMINACLIP:${params.trimmomatic_adater_path}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    """
}

process nanofilt {
    label "containerPython"
    conda "/home/molecularvirology/miniconda2/envs/vdp_lrs"
    input:
        path fastq
    output:
        path "${params.prefix}.filtered.fastq"
    script:
        
    """
    echo "Filter low quality reads"
    echo "sample name: ${params.prefix}"
    echo "input file: $fastq"
    
    NanoFilt --length 300 --quality 9 --readtype 1D $fastq > ${params.prefix}.filtered.fastq

    """
}

process hostfilter_illumina {
//    errorStrategy { 'ignore' }

    label "containerBowtie2"
    conda "/home/molecularvirology/miniconda2/envs/vdp_srs"

    input:
        tuple path(pe1), path(pe2)
    output:
        tuple path("${params.prefix}_hostfiltered_1.fastq"), path("${params.prefix}_hostfiltered_2.fastq")
    """
    echo "Host filter with minimap2 and samtools"
    echo "input file 1: $pe1"
    echo "input file 2: $pe2"
    echo "host organism: ${params.host}"
    bowtie2 -x ${params.host_ref_path} -1 $pe1 -2 $pe2 -S temp.sam
    samtools view -Sb temp.sam | samtools view -b -f 12 -F 256 | samtools sort -o  "${params.prefix}.hostfiltered.bam"
    rm temp.sam
    bedtools bamtofastq -i "${params.prefix}.hostfiltered.bam" -fq ${params.prefix}_hostfiltered_1.fastq -fq2 ${params.prefix}_hostfiltered_2.fastq
    """
}

process hostfilter_nanopore {
//    errorStrategy { 'ignore' }

    label "containerHostfilter"
    conda "/home/molecularvirology/miniconda2/envs/vdp_lrs"
    publishDir "${params.outdir}/filter", mode: 'copy'
    
    input:
        path fastq
    output:
        path "${params.prefix}.hostfiltered.fastq"
    """
    echo "Host filter with minimap2 and samtools"
    echo "input file: $fastq"
    echo "host organism: $params.host"
    minimap2 -ax map-ont $params.host_ref_path $fastq | samtools view -Sb - | samtools sort -o aln.sorted.bam
    samtools fastq -f 4 aln.sorted.bam > "${params.prefix}.hostfiltered.fastq"
    """
}
