nextflow.enable.dsl=2

workflow classify_reads_illumina {
    take:
        filtered_fastqs
    main:
        kraken_result = kraken2_illumina(filtered_fastqs)
        kreport2krona(kraken_result)
        edge_form = kreport2EDGEform(kraken_result)
        EDGEform2heatmap(edge_form)
        //kaiju_results = kaiju_illumina(filtered)
        //kaiju2krona(kaiju_results)
}

workflow classify_reads_nanopore {
    take:
        filtered
    main:
        kraken_result = kraken2_nanopore(filtered)
        kreport2krona(kraken_result)
        edge_form = kreport2EDGEform(kraken_result)
        EDGEform2heatmap(edge_form)
        kaiju_result = kaiju_nanopore(filtered)
        kaiju2krona(kaiju_result)
}

process kreport2EDGEform {
    conda "/home/molecularvirology/miniconda2/envs/vdp_lrs"
    label "containerKrona"

    input:
        path kraken_report
    output:
        path "${params.prefix}.list"

    """
    perl ${params.pipeline_directory}/scripts/convert_krakenRep2list.pl < $kraken_report > ${params.prefix}.list
    """
}

process EDGEform2heatmap {
    publishDir "${params.outdir}/classification", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/vdp_srs"
    label "containerMetaComp"

    input:
        path edge_form
    output:
        path "${params.prefix}_phylum.svg"
        path "${params.prefix}_family.svg"
        path "${params.prefix}_species.svg"
    """
    Rscript ${params.pipeline_directory}/scripts/make_heatmap_with_metacomp.R ${edge_form} ${params.prefix}
    """
}

process kraken2_illumina {
    conda "/home/molecularvirology/miniconda2/envs/vdp_lrs"
    label "containerKraken"
    input:
        tuple path(pe1), path(pe2)
    output:
        path "${params.prefix}.kraken_report.csv"
    """
    echo "Classify reads into taxonomy with Kraken2"
    echo "sample name: $params.prefix"
    echo "input file 1: $pe1"
    echo "input file 2: $pe2"
    echo "kraken db path: $params.kraken_db_path"
    kraken2 --db $params.kraken_db_path \
        --report ${params.prefix}.kraken_report.csv \
        --paired --threads 12 \
        $pe1 $pe2
    """
}

process kraken2_nanopore {
    conda "/home/molecularvirology/miniconda2/envs/vdp_lrs"
    label "containerKraken"
    input:
        path fastq
    output:
        path "${params.prefix}.kraken_report.csv"
    """
    echo "Classify reads into taxonomy with Kraken2"
    echo "sample name: $params.prefix"
    echo "query file: $fastq"
    echo "kraken db path: $params.kraken_db_path"
    kraken2 --db $params.kraken_db_path \
        --report ${params.prefix}.kraken_report.csv \
        --threads 12 $fastq
    """
}

process kreport2krona {
    publishDir "${params.outdir}/classification", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/vdp_lrs"
    label "containerKrona"
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
    publishDir "${params.outdir}/classification", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/vdp_lrs"
    label "containerKaiju"
    input:
        path fastq
    output:
        path "${params.prefix}.kaiju.out.krona"
    """
    echo "Classify reads into taxonomy with Kaiju"
    echo "sample name: $params.prefix"
    echo "query file: $fastq"
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
    publishDir "${params.outdir}/classification", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/vdp_lrs"
    label "containerKrona"
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

process bracken {
    errorStrategy 'ignore'
    publishDir "${params.outdir}/classification", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/vdp_lrs"
    label "containerBracken"
    input:
        path kraken_report
    output:
        path "${params.prefix}.bracken_all.txt"

    script:
    """
    echo "Estimate abundance of taxons at a single level(default species) with bracken"
    echo "sample name: $params.prefix"

    bracken -d $params.kraken_db_path -i $kraken_report -o ${params.prefix}.bracken_all.txt -l 'S' -t 10 -r 200
    """
}
