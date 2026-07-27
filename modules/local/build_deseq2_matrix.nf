process BUILD_DESEQ2_MATRIX {
    tag "$tier"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/21497240/files/nanocirc-quant-latest.sif?download=1' :
        'quay.io/anrusakovich/nanocirc-quant:latest' }"

    input:
    tuple val(tier), val(metas), path(counts_tsvs)

    output:
    tuple val(tier), path("deseq2_counts_${tier}.tsv"),   emit: counts
    tuple val(tier), path("deseq2_coldata_${tier}.tsv"),  emit: coldata
    tuple val(tier), path("deseq2_features_${tier}.tsv"), emit: features
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample_args = [metas, counts_tsvs].transpose().collect { meta, tsv -> "${meta.id}=${tsv}" }.join(' ')
    def group_args  = metas.collect { meta -> "${meta.id}=${meta.group ?: meta.id}" }.join(' ')
    """
    build_deseq2_matrix.py \\
        --tier ${tier} \\
        --counts_tsv ${sample_args} \\
        --groups ${group_args} \\
        --out_counts   deseq2_counts_${tier}.tsv \\
        --out_coldata  deseq2_coldata_${tier}.tsv \\
        --out_features deseq2_features_${tier}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
