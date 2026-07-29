process QUANT_APPEND_COUNTS {
    tag "${meta.id}:${meta.category}"
    label 'process_single'
    container "${ (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/21497240/files/nanocirc-quant-latest.sif?download=1' :
        'quay.io/anrusakovich/nanocirc-quant:latest' }"

    input:
    tuple val(meta), path(clean_tsv), path(cluster_map), path(quant_counts)

    output:
    tuple val(meta), path("${meta.id}_${meta.category}_clean_with_counts.tsv"), emit: clean
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quant_append_counts.py \\
        --clean_tsv ${clean_tsv} \\
        --cluster_map ${cluster_map} \\
        --quant_counts ${quant_counts} \\
        --output ${meta.id}_${meta.category}_clean_with_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
