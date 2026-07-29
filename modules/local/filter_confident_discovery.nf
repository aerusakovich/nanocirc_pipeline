process FILTER_CONFIDENT_DISCOVERY {
    tag "${meta.id}:${meta.category}"
    label 'process_single'
    container "${ (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/21497240/files/nanocirc-quant-latest.sif?download=1' :
        'quay.io/anrusakovich/nanocirc-quant:latest' }"

    stageInMode 'copy'

    input:
    tuple val(meta), path(clean_with_counts_tsv), path(unfiltered_bed)

    output:
    tuple val(meta), path("${meta.id}_${meta.category}_clean_with_counts.tsv"), emit: clean
    tuple val(meta), path("${meta.id}_${meta.category}.bed12"),                 emit: bed
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_confident_discovery.py \\
        --tsv       ${clean_with_counts_tsv} \\
        --bed       ${unfiltered_bed} \\
        --min_reads ${params.circrna_confident_min_reads} \\
        --category  ${meta.category} \\
        --out_tsv   ${meta.id}_${meta.category}_clean_with_counts.tsv \\
        --out_bed   ${meta.id}_${meta.category}.bed12

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
