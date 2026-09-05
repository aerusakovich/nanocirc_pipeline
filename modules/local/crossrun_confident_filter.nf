process CROSSRUN_CONFIDENT_FILTER {
    tag "${meta.id}:${meta.tier}"
    label 'process_single'
    container "${ (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/21497240/files/nanocirc-quant-latest.sif?download=1' :
        'quay.io/anrusakovich/nanocirc-quant@sha256:11534ee495f3717000a9b848600a395db092d22798a40eb4bccd533f7ca362ad' }"

    stageInMode 'copy'

    input:
    tuple val(meta), path(group_clean_tsv), path(group_bed12), path(run_filtered_tsvs)

    output:
    tuple val(meta), path("${prefix}_clean.tsv"), emit: clean
    tuple val(meta), path("${prefix}.bed12"),      emit: bed
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${meta.id}_${meta.tier}_crossrun"
    """
    crossrun_confident_filter.py \\
        --group_clean_tsv    ${group_clean_tsv} \\
        --group_bed12        ${group_bed12} \\
        --run_filtered_tsvs  ${run_filtered_tsvs} \\
        --out_tsv             ${prefix}_clean.tsv \\
        --out_bed12           ${prefix}.bed12

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
