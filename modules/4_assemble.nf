nextflow.enable.dsl=2

workflow assemble_illumina {
    take:
        fastq
    emit:
        renamed_contigs
    main:
        contigs = spades(fastq)
        filtered_contigs = filterContigs_illumina(contigs)
        renamed_contigs = renameContigs(filtered_contigs)
        contigLen(renamed_contigs)
        
}

workflow assemble_nanopore {
    take:
        fastq
    emit:
        renamed_contigs
    main:
        contigs = megahit(fastq)
        filtered_contigs = filterContigs_nanopore(contigs)
        renamed_contigs = renameContigs(filtered_contigs)
        contigLen(renamed_contigs)
}

process filterContigs_illumina {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'

    input:
        path contigs
    output:
        path("${params.prefix}.filtered_contigs.fasta")
    """
    reformat.sh in=$contigs out="${params.prefix}.filtered_contigs.fasta" minlength=${params.illumina_min_contig_length}
    """
}

process filterContigs_nanopore {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    input:
        path contigs
    output:
        path("${params.prefix}.filtered_contigs.fasta")
    """
    reformat.sh in=$contigs out="${params.prefix}.filtered_contigs.fasta" minlength=${params.nanopore_min_contig_length}
    """
}

process renameContigs {
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
        path contigs
    output:
        path("${params.prefix}.contigs.fasta")
    """
    awk '/>/{print ">tig" ++i; next}{print}' < $contigs > "${params.prefix}.contigs.fasta"
    """
}

process contigLen {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:
        path contigs
    output:
        path "${params.prefix}.contigs.len"
    """
    python ${params.nextflow_script_path}/scripts/4_contigs_length.py $params.prefix $contigs
    """
}

process spades{
    conda '/home/molecularvirology/miniconda2/envs/vdp_srs'
    errorStrategy { 'ignore' }
    publishDir "${params.outdir}/assembly", mode: 'copy'
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
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
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
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    errorStrategy { 'ignore' }
    publishDir "${params.outdir}/assembly", mode: 'copy'
    
    input:
        path fastq
    output:
        path("${params.prefix}/${params.prefix}.contigs.fasta") optional true
    """
    echo "De novo assembly with Canu"
    canu -p $params.prefix -d $params.prefix -nanopore $fastq \
    genomeSize=5m minReadLength=300 minOverlapLength=50 maxThreads=12 \
    minInputCoverage=0 stopOnLowCoverage=0 corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 \
    redMemory=32 oeaMemory=32 batMemory=200 useGrid=false
    """
}
