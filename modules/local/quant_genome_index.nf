process QUANT_GENOME_INDEX {
    tag "genome_index"
    label 'process_high'
    errorStrategy 'retry'
    maxRetries    2

    // Cache results, skipped on subsequent runs if the index already exists
    // here, same pattern as PREPARE_GENOME.
    // retry: storeDir unstage/verify race (nxf_unstage moves output before
    // Nextflow's own existence check reads it)
    storeDir "${params.genome_index_dir ?: "${params.outdir}/genome_index"}"

    container "${ (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/21497240/files/nanocirc-quant-latest.sif?download=1' :
        'quay.io/anrusakovich/nanocirc-quant@sha256:11534ee495f3717000a9b848600a395db092d22798a40eb4bccd533f7ca362ad' }"

    input:
    path fasta

    output:
    path "genome_splice.mmi", emit: index
    path "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quant_genome_index.py --genome_fasta ${fasta} --output genome_splice.mmi --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """
}
