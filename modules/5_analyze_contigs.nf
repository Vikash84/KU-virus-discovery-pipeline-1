nextflow.enable.dsl=2

workflow analyze_contigs {
    take:
        contigs
    main:
        blastn_results = blastn(contigs)
        blastn_aln = blastn_results[0]
        megablast_aln = blastn_results[1]

        makeHtmlFromBlastAln_1(blastn_aln)
        blastn_tmp = makeTxtFromBlastAln_1(blastn_aln)
        blastn_txt = add_header_1(blastn_tmp)
        filterBlastResult_1(blastn_txt)

        makeHtmlFromBlastAln_2(megablast_aln)
        megablast_tmp = makeTxtFromBlastAln_2(megablast_aln)
        megablast_txt = add_header_2(megablast_tmp)
        filterBlastResult_2(megablast_txt)

        blastx_tmp = diamondBlastx(contigs)
        blastx_txt = add_header_3(blastx_tmp)
        filterBlastResult_3(blastx_txt)
}

process blastn {
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerBlast"
    
    input:
        path contigs
    output:
        path "${params.prefix}.blastn.aln"
        path "${params.prefix}.megablast.aln"

    """
    blastn -query $contigs -db $params.blastn_db_path -task blastn -evalue 1.0e-5 -max_target_seqs 1 -outfmt 11 > ${params.prefix}.blastn.aln
    blastn -query $contigs -db $params.blastn_db_path -task megablast -evalue 1.0e-5 -max_target_seqs 1 -outfmt 11 > ${params.prefix}.megablast.aln
   """
}

process diamondBlastx {
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerDiamond"

    input:
        path contigs
    output:
        path "${params.prefix}.blastx.tmp"
    """
    diamond blastx -d $params.diamond_db_path -q $contigs --max-target-seqs 1 --evalue 1.0e-5 --outfmt 6 qseqid sseqid stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore > ${params.prefix}.blastx.tmp
    """
}

process makeHtmlFromBlastAln_1 {
    publishDir "${params.outdir}/analysis", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerBlast"

    input:
        path aln
    output:
        path "${aln.baseName}.html"
    """
    blast_formatter -archive $aln -html > ${aln.baseName}.html
    """
}

process makeHtmlFromBlastAln_2 {
    publishDir "${params.outdir}/analysis", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerBlast"

    input:
        path aln
    output:
        path "${aln.baseName}.html"
    """
    blast_formatter -archive $aln -html > ${aln.baseName}.html
    """
}

process makeTxtFromBlastAln_1 {
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerBlast"

    input:
        path aln
    output:
        path "${aln.baseName}.tmp"
    """
    blast_formatter -archive $aln -outfmt "6 qseqid sseqid stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore"  > "${aln.baseName}.tmp"
    """
}

process makeTxtFromBlastAln_2 {
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerBlast"

    input:
        path aln
    output:
        path "${aln.baseName}.tmp"
    """
    blast_formatter -archive $aln -outfmt "6 qseqid sseqid stitle pident qlen length mismatch gapopen qstart qend sstart send evalue bitscore"  > "${aln.baseName}.tmp"
    """
}

process add_header_1 {
    input:
        path tmp
    output:
        path "${tmp.baseName}.txt"
    """
    cat ${params.pipeline_directory}/headers/blast_header $tmp > ${tmp.baseName}.txt
    """
}

process add_header_2 {
    input:
        path tmp
    output:
        path "${tmp.baseName}.txt"
    """
    cat ${params.pipeline_directory}/headers/blast_header $tmp > ${tmp.baseName}.txt
    """
}

process add_header_3 {
    input:
        path tmp
    output:
        path "${tmp.baseName}.txt"
    """
    cat ${params.pipeline_directory}/headers/blast_header $tmp > ${tmp.baseName}.txt
    """
}

process filterBlastResult_1 {
    publishDir "${params.outdir}/analysis", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerPython"

    input:
        path txt
    output:
        path "${txt.baseName}.filtered.txt"

    """
    python ${params.pipeline_directory}/scripts/blast_filter.py --input $txt --output ${txt.baseName}.filtered.txt --exclude ${params.exclude_keyword_file} --aln_len 100 --aln_contig_len_prop 0.5
    """
}

process filterBlastResult_2 {
    publishDir "${params.outdir}/analysis", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerPython"

    input:
        path txt
    output:
        path "${txt.baseName}.filtered.txt"

    """
    python ${params.pipeline_directory}/scripts/blast_filter.py --input $txt --output ${txt.baseName}.filtered.txt --exclude ${params.exclude_keyword_file} --aln_len 100 --aln_contig_len_prop 0.5
    """
}

process filterBlastResult_3 {
    publishDir "${params.outdir}/analysis", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerPython"

    input:
        path txt
    output:
        path "${txt.baseName}.filtered.txt"

    """
    python ${params.pipeline_directory}/scripts/blast_filter.py --input $txt --output ${txt.baseName}.filtered.txt --exclude ${params.exclude_keyword_file} --aln_len 33 --aln_contig_len_prop 0.17
    """
}

process vibrant {
    publishDir "${params.outdir}/analysis", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerVibrant"

    input:
        path contigs
    output:
        path "VIBRANT/**VIBRANT_summary_results_*" optional true
    """
    VIBRANT_run.py -i $contigs -folder VIBRANT -d ${params.vibrant_db_path} -m ${params.vibrant_files_path}
    """
}

process deepvirfinder {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/analysis", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerDVF"

    input:
        path contigs
    output:
        path "DVF/*dvfpred.txt"
    """
    python ${params.deepvirfinder_path}/dvf.py -i $contigs -o DVF -l 1000 -c 16
    Rscript ${params.pipeline_directory}/scripts/dvf_filter.R DVF/*dvfpred.txt
    """
}

process parse_and_prokka {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/analysis", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerBlast"

    input:
        path dvf_pred
        path contigs
    output:
        path "prokka/${params.prefix}.tsv"
        path "prokka/${params.prefix}.gff"
    """
    python ${params.pipeline_directory}/scripts/parse_dvf_pred.py $dvf_pred $contigs ${params.prefix}.contigs
    prokka --kingdom Viruses --outdir prokka --params.prefix ${params.prefix} ${params.prefix}.contigs
    """
}
