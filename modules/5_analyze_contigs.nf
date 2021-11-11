nextflow.enable.dsl=2

workflow analyze_contigs {
    take:
        contigs
    main:
        blastn_results = blastn(contigs)
        blastn_result = blastn_results[0]
        megablast_result = blastn_results[1]
        diamond_result = diamondBlastx(contigs)

        parseBlastn(blastn_result)
        parseMegaBlast(megablast_result)
        parseDiamondBlastx(diamond_result)
}

process blastn {
    
    input:
        path contigs
    output:
        path "${params.prefix}.blastn.tmp"
        path "${params.prefix}.megablast.tmp"

    """
    blastn -query $contigs -db $params.blastn_db_path -task blastn -evalue 1.0e-5 -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 12 > ${params.prefix}.blastn.tmp
    blastn -query $contigs -db $params.blastn_db_path -task megablast -evalue 1.0e-5 -max_target_seqs 1 -outfmt "6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 12 > ${params.prefix}.megablast.tmp
   """
}

process diamondBlastx {

    input:
        path contigs
    output:
        path "${params.prefix}.blastx.tmp"
    """
    diamond blastx -d $params.diamond_db_path -q $contigs --max-target-seqs 1 --evalue 1.0e-5 --outfmt 6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore > ${params.prefix}.blastx.tmp
    """
}

process parseBlastn {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/analysis", mode: 'copy'

    input:
        path result
    output:
        path "${params.prefix}.blastn.txt"
    """
    cat ${params.nextflow_script_path}/headers/blast_header $result > ${result.baseName}.headered.tmp
    python ${params.nextflow_script_path}/scripts/5_blast_result_filter.py --input ${result.baseName}.headered.tmp --output ${result.baseName}.filtered.tmp --aln_len 100 --aln_contig_len_prop 0.5
    Rscript ${params.nextflow_script_path}/scripts/5_match_taxonomy_to_blast_result.R ${result.baseName}.filtered.tmp ${params.prefix}.blastn.txt ${params.taxonomizr_db_path}
    """
}

process parseMegaBlast {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/analysis", mode: 'copy'

    input:
        path result
    output:
        path "${params.prefix}.megablast.txt"
    """
    cat ${params.nextflow_script_path}/headers/blast_header $result > ${result.baseName}.headered.tmp
    python ${params.nextflow_script_path}/scripts/5_blast_result_filter.py --input ${result.baseName}.headered.tmp --output ${result.baseName}.filtered.tmp --aln_len 100 --aln_contig_len_prop 0.5
    Rscript ${params.nextflow_script_path}/scripts/5_match_taxonomy_to_blast_result.R ${result.baseName}.filtered.tmp ${params.prefix}.megablast.txt ${params.taxonomizr_db_path}
    """
}

process parseDiamondBlastx {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/analysis", mode: 'copy'

    input:
        path result
    output:
        path "${params.prefix}.blastx.txt"
    """
    cat ${params.nextflow_script_path}/headers/blast_header $result > ${result.baseName}.headered.tmp
    python ${params.nextflow_script_path}/scripts/5_blast_result_filter.py --input ${result.baseName}.headered.tmp --output ${result.baseName}.filtered.tmp --aln_len 33 --aln_contig_len_prop 0.17
    Rscript ${params.nextflow_script_path}/scripts/5_match_taxonomy_to_blast_result.R ${result.baseName}.filtered.tmp ${params.prefix}.blastx.txt ${params.taxonomizr_db_path}
    """
}