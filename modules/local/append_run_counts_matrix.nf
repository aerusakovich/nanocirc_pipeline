process APPEND_RUN_COUNTS_MATRIX {
    tag "$meta.id"
    label 'process_single'
    container "${ (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/21497240/files/nanocirc-quant-latest.sif?download=1' :
        'quay.io/anrusakovich/nanocirc-quant@sha256:11534ee495f3717000a9b848600a395db092d22798a40eb4bccd533f7ca362ad' }"

    stageInMode 'copy'

    input:
    tuple val(meta), path(clean_tsv), path(counts_matrix)

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: clean
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.category}_clean_with_counts"
    """
    append_run_counts_matrix.py \\
        --clean_tsv ${clean_tsv} \\
        --counts_matrix ${counts_matrix} \\
        --output ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
