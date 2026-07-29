process QUANT_COVERAGE_FLAG {
    tag "$meta.id"
    label 'process_medium'
    container "${ (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/21497240/files/nanocirc-quant-latest.sif?download=1' :
        'quay.io/anrusakovich/nanocirc-quant:latest' }"

    input:
    tuple val(meta), path(sample_bed12), path(sample_conf_tsv),
          path(isocirc_native), path(circfl_native), path(cirilong_reads), path(circnick_native),
          path(unit_deduped_metadata), path(unit_cluster_map), path(tier1_counts)

    output:
    tuple val(meta), path("${meta.id}_flagged_loci.tsv"), emit: flagged
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_old_tool_count = params.quant_min_old_tool_count ?: 50
    // Inactive or absent tool inputs arrive as a per-tool NO_FILE_* sentinel
    // (see circrna_analysis.nf). Only pass the corresponding flag through
    // when a real file was staged.
    def iso_arg  = isocirc_native.name  != 'NO_FILE_ISOCIRC'  ? "--isocirc_out ${isocirc_native}"        : ''
    def fl_arg   = circfl_native.name   != 'NO_FILE_CIRCFL'   ? "--circfl_pass ${circfl_native}"         : ''
    def ciri_arg = cirilong_reads.name  != 'NO_FILE_CIRILONG' ? "--cirilong_reads ${cirilong_reads}"     : ''
    def nick_arg = circnick_native.name != 'NO_FILE_CIRCNICK' ? "--circnick_annotated ${circnick_native}" : ''
    """
    quant_catalog_from_bed12.py --bed ${sample_bed12} --tsv ${sample_conf_tsv} --output ${meta.id}_sample_catalog.tsv

    quant_tool_readset_union.py \\
        --catalog ${meta.id}_sample_catalog.tsv \\
        --tsv ${sample_conf_tsv} \\
        --sample ${meta.id} \\
        ${iso_arg} ${fl_arg} ${ciri_arg} ${nick_arg}

    quant_coverage_flag.py \\
        --deduped_metadata ${unit_deduped_metadata} \\
        --cluster_map ${unit_cluster_map} \\
        --remap_counts ${tier1_counts} \\
        --old_tool_counts ${meta.id}_old_tool_counts.tsv \\
        --sample ${meta.id} \\
        --min_old_tool_count ${min_old_tool_count}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
