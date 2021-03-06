nextflow.enable.dsl=2
include { qc_nanopore } from './modules/0_qc'
include { filter_nanopore } from './modules/1_filter'
include { reference_mapping_nanopore } from './modules/2_reference_mapping'
include { classify_reads_nanopore } from './modules/3_classify_reads'
include { assemble_nanopore } from './modules/4_assemble'
include { polish } from './modules/4_2_polish'
include { analyze_contigs } from './modules/5_analyze_contigs'

if(params.help) {
    log.info ''
    log.info 'KU virus discovery pipeline - long read (Nanopore) version'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run main_nanopore.nf [options]'
    log.info ''
    log.info 'Script Options: '
    log.info '    --fastq               FILE    Path to FASTQ file'
    log.info '    --prefix              STR     Nickname given to sample'
    log.info '    --outdir              DIR     Name of output directory'
    log.info '    --host                STR     Name of host the sample originate from'
    log.info '    --host_reference_path STR     Path to the fast of host genome'
    log.info ''

    return
}

// workflow module
workflow {
    main:
        fastq = channel.fromPath(params.fastq, checkIfExists:true)

        // fastq quality control
        qc_nanopore(fastq)

        // filter too short/low quality and host-derived reads
        filtered = filter_nanopore(fastq)

        // reference mapping with provided virus sequence lists
        reference_mapping_nanopore(filtered)

        // read taxon classification
        classify_reads_nanopore(filtered)

        // de novo assembly
        contigs = assemble_nanopore(filtered)

        // post-assembly work
        polished_contigs = polish(contigs, fastq)

        // contigs homology search and functional analysis
        analyze_contigs(polished_contigs)
}
