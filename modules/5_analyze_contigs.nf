nextflow.enable.dsl=2

workflow analyze_contigs {
    take:
        contigs
    main:
        blastn_results = blastn(contigs)
        blastn_aln = blastn_results[0]
        megablast_aln = blastn_results[1]

        blastn_tmp = makeTxtFromBlastAln_1(blastn_aln)
        blastn_txt = add_header_1(blastn_tmp)
        blastn_filtered = filterBlastResult_blastn(blastn_txt)
        matchTaxonomyToBlastResult_blastn(blastn_filtered)

        megablast_tmp = makeTxtFromBlastAln_2(megablast_aln)
        megablast_txt = add_header_2(megablast_tmp)
        megablast_filtered = filterBlastResult_megablast(megablast_txt)
        matchTaxonomyToBlastResult_megablast(megablast_filtered)

        blastx_tmp = diamondBlastx(contigs)
        blastx_txt = add_header_3(blastx_tmp)
        blastx_filtered = filterBlastResult_blastx(blastx_txt)
        matchTaxonomyToBlastResult_blastx(blastx_filtered)
}

process blastn {
    
    input:
        path contigs
    output:
        path "${params.prefix}.blastn.aln"
        path "${params.prefix}.megablast.aln"

    """
    blastn -query $contigs -db $params.blastn_db_path -task blastn -evalue 1.0e-5 -max_target_seqs 1 -outfmt 11 -num_threads 12 > ${params.prefix}.blastn.aln
    blastn -query $contigs -db $params.blastn_db_path -task megablast -evalue 1.0e-5 -max_target_seqs 1 -outfmt 11 -num_threads 12 > ${params.prefix}.megablast.aln
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

process makeTxtFromBlastAln_1 {

    input:
        path aln
    output:
        path "${aln.baseName}.tmp"
    """
    blast_formatter -archive $aln -outfmt "6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore"  > "${aln.baseName}.tmp"
    """
}

process makeTxtFromBlastAln_2 {

    input:
        path aln
    output:
        path "${aln.baseName}.tmp"
    """
    blast_formatter -archive $aln -outfmt "6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore"  > "${aln.baseName}.tmp"
    """
}

process add_header_1 {
    input:
        path tmp
    output:
        path "${tmp.baseName}.headered.tmp"
    """
    cat ~/headers/blast_header $tmp > ${tmp.baseName}.headered.tmp
    """
}

process add_header_2 {
    input:
        path tmp
    output:
        path "${tmp.baseName}.headered.tmp"
    """
    cat ~/headers/blast_header $tmp > ${tmp.baseName}.headered.tmp
    """
}

process add_header_3 {
    input:
        path tmp
    output:
        path "${tmp.baseName}.headered.tmp"
    """
    cat ~/headers/blast_header $tmp > ${tmp.baseName}.headered.tmp
    """
}

process filterBlastResult_blastn {

    input:
        path txt
    output:
        path "${txt.baseName}.filtered.tmp"

    """
    python ~/scripts/5_blast_result_filter.py --input $txt --output ${txt.baseName}.filtered.tmp --aln_len 100 --aln_contig_len_prop 0.5
    """
}

process filterBlastResult_megablast {

    input:
        path txt
    output:
        path "${txt.baseName}.filtered.tmp"

    """
    python ~/scripts/5_blast_result_filter.py --input $txt --output ${txt.baseName}.filtered.tmp --aln_len 100 --aln_contig_len_prop 0.5
    """
}

process filterBlastResult_blastx {

    input:
        path txt
    output:
        path "${txt.baseName}.filtered.tmp"

    """
    python ~/scripts/5_blast_result_filter.py --input $txt --output ${txt.baseName}.filtered.tmp --aln_len 33 --aln_contig_len_prop 0.17
    """
}

process matchTaxonomyToBlastResult_blastn {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/analysis", mode: 'copy'

    input:
        path txt
    output:
        path "${txt.simpleName}.blastn.txt" optional true

    """
    Rscript ~/scripts/5_match_taxonomy_to_blast_result.R $txt ${txt.simpleName}.blastn.txt ${params.taxonomizr_db_path}
    """
}

process matchTaxonomyToBlastResult_megablast {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/analysis", mode: 'copy'

    input:
        path txt
    output:
        path "${txt.simpleName}.megablast.txt" optional true

    """
    Rscript ~/scripts/match_taxonomy_to_blast.R $txt ${txt.simpleName}.megablast.txt ${params.taxonomizr_db_path}
    """
}

process matchTaxonomyToBlastResult_blastx {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/analysis", mode: 'copy'

    input:
        path txt
    output:
        path "${txt.simpleName}.blastx.txt" optional true

    """
    Rscript ~/scripts/match_taxonomy_to_blast.R $txt ${txt.simpleName}.blastx.txt ${params.taxonomizr_db_path}
    """
}

