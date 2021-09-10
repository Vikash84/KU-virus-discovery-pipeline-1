nextflow.enable.dsl=2

workflow classify_reads {
    take:
        filtered
    main:
        kraken_results = kraken2(filtered)
        kreport2krona(kraken_results[0])
        kaiju_results = kaiju(filtered)
        kaiju2krona(kaiju_results[2])
}

process kraken2 {
    publishDir "${params.outDir}/classification", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerKraken"
    input:
        path fastq
    output:
        path "${params.prefix}.kraken_report.txt"
        path "${params.prefix}.kraken_output.txt"
    """
    echo "Classify reads into taxonomy with Kraken2"
    echo "sample name: $params.prefix"
    echo "query file: $fastq"
    echo "kraken db path: $params.kraken_db_path"
    kraken2 --db $params.kraken_db_path \
        --report ${params.prefix}.kraken_report.txt \
        --output ${params.prefix}.kraken_output.txt \
        $fastq
    """
}

process kreport2krona {
    publishDir "${params.outDir}/classification", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
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

process kaiju {
    publishDir "${params.outDir}/classification", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
    label "containerKaiju"
    input:
        path fastq
    output:
        path "${params.prefix}.kaiju.out"
        path "${params.prefix}.kaiju_summary.tsv"
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
    publishDir "${params.outDir}/classification", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
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
    publishDir "${params.outDir}/classification", mode: 'copy'
    conda "/home/molecularvirology/miniconda2/envs/ku_vdp"
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
