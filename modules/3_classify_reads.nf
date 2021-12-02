nextflow.enable.dsl=2

workflow classify_reads_illumina {
    take:
        fastqs
    main:
        kraken_result = kraken2_illumina(fastqs)
        kreport2krona(kraken_result)
        edge_form = kreport2EDGEform(kraken_result)
        EDGEform2heatmap(edge_form)
}

workflow classify_reads_nanopore {
    take:
        fastq
    main:
        kraken_result = kraken2_nanopore(fastq)
        kreport2krona(kraken_result)
        edge_form = kreport2EDGEform(kraken_result)
        EDGEform2heatmap(edge_form)

        kaiju_result = kaiju_nanopore(fastq)
        kaiju2krona(kaiju_result)
}

process kreport2EDGEform {
    conda '/home/molecularvirology/miniconda2/envs/vdp_srs'
    input:
        path kraken_report
    output:
        path "${params.prefix}.list"

    """
    perl ${params.nextflow_script_path}/scripts/3_convert_krakenRep2list.pl < $kraken_report > ${params.prefix}.list
    """
}

process EDGEform2heatmap {
    conda '/home/molecularvirology/miniconda2/envs/vdp_srs'
    publishDir "${params.outdir}/classification", mode: 'copy'

    input:
        path edge_form
    output:
        path "${params.prefix}_order.svg"
        path "${params.prefix}_family.svg"
        path "${params.prefix}_genus.svg"
        path "${params.prefix}_species.svg"
    """
    Rscript ${params.nextflow_script_path}/scripts/3_generate_heatmap.R ${edge_form} ${params.prefix}
    """
}

process kraken2_illumina {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    input:
        tuple path(pe1), path(pe2)
    output:
        path "${params.prefix}.kraken_report.csv"
    """
    echo "Classify reads into taxonomy with Kraken2"
    echo "sample name: $params.prefix"
    echo "kraken db path: $params.kraken_db_path"
    kraken2 --db $params.kraken_db_path \
        --report ${params.prefix}.kraken_report.csv \
        --paired --threads 12 --confidence 0.1\
        $pe1 $pe2
    """
}

process kraken2_nanopore {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    input:
        path fastq
    output:
        path "${params.prefix}.kraken_report.csv"
    """
    echo "Classify reads into taxonomy with Kraken2"
    echo "sample name: $params.prefix"
    echo "kraken db path: $params.kraken_db_path"
    kraken2 --db $params.kraken_db_path \
        --report ${params.prefix}.kraken_report.csv \
        --threads 12 --confidence 0.1\
        $fastq
    """
}

process kreport2krona {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    publishDir "${params.outdir}/classification", mode: 'copy'
    input:
        path kraken_report
    output:
        path "${params.prefix}.kraken.html"
    """
    echo "Convert Kraken report to Krona format"
    echo "sample name: $params.prefix"
    kreport2krona.py -r $kraken_report -o ${params.prefix}.krona
    ktImportText ${params.prefix}.krona -o ${params.prefix}.kraken.html
    """
}

process kaiju_nanopore {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    publishDir "${params.outdir}/classification", mode: 'copy'
    input:
        path fastq
    output:
        path "${params.prefix}.kaiju.out.krona"
    """
    echo "Classify reads into taxonomy with Kaiju"
    echo "sample name: $params.prefix"
    echo "kaiju db path: $params.kaiju_db_path"
    kaiju -f $params.kaiju_db_path/kaiju_db_nr.fmi \
        -t $params.kaiju_db_path/nodes.dmp \
        -i $fastq \
        -o ${params.prefix}.kaiju.out \
        -a greedy -E 0.05 -z 24 -v
    kaiju2table -t ${params.kaiju_db_path}/nodes.dmp -n ${params.kaiju_db_path}/names.dmp -r species -e -o ${params.prefix}.kaiju_summary.tsv ${params.prefix}.kaiju.out
    kaiju2krona -t ${params.kaiju_db_path}/nodes.dmp -n ${params.kaiju_db_path}/names.dmp -i ${params.prefix}.kaiju.out -o ${params.prefix}.kaiju.out.krona
    """
}

process kaiju2krona {
    conda '/home/molecularvirology/miniconda2/envs/vdp_lrs'
    publishDir "${params.outdir}/classification", mode: 'copy'
    input:
        path kaiju_krona
    output:
        path "${params.prefix}.kaiju.html"
    """
    echo "Visualize Kaiju result with Krona"
    echo "sample name: $params.prefix"
    ktImportText $kaiju_krona -o ${params.prefix}.kaiju.html
    """
}