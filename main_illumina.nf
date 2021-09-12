nextflow.enable.dsl=2
include { qc_illumina } from './modules/0_qc'
include { filter_illumina } from './modules/1_filter'
include { reference_map_illumina } from './modules/2_reference_map'
include { classify_reads_illumina } from './modules/3_classify_reads'
include { assemble_illumina } from './modules/4_assemble'
include { analyze_contigs } from './modules/5_analyze_contigs'

if(params.help) {
    log.info ''
    log.info 'KU virus discovery pipeline - short read (Illumina) version'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run main.nf [options]'
    log.info ''
    log.info 'Script Options: '
    log.info '    --fastq        FILE    Path to first paired-end FASTQ file'
    log.info '    --fastq2        FILE    Path to second paired-end FASTQ file'
    log.info '    --prefix        STR    Nickname given to sample'
    log.info '    --outdir         DIR      Name of output directory'
    log.info ''

    return
}


// workflow module
workflow {
    main:
        fastq1 = channel.fromPath(params.fastq, checkIfExists:true)
        fastq2 = channel.fromPath(params.fastq2, checkIfExists:true)

        // fastq quality control
//        qc_illumina(fastq1, fastq2)

        // filter too short/low quality and host-derived reads
        filtered = filter_illumina(fastq1, fastq2)

        // reference mapping with provided virus sequence lists
//        reference_map_illumina(filtered)

        // read taxon classification
//        classify_reads_illumina(filtered)

        // de novo assembly
        contigs = assemble_illumina(filtered)

        // contigs homology search and functional analysis
        analyze_contigs(contigs)
}
