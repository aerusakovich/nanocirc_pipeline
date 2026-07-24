process QUANT_PREPARE_CATALOG {
    tag "$unit_id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/21497240/files/nanocirc-quant-latest.sif?download=1' :
        'quay.io/anrusakovich/nanocirc-quant:latest' }"

    input:
    tuple val(unit_id), path(bed12), path(conf_tsv)
    path fasta

    output:
    tuple val(unit_id), path("${unit_id}_deduped_metadata.tsv"), path("${unit_id}_cluster_map.tsv"),
          path("circle_refs.fa"), path("circle_refs_lengths.tsv"), path("chunks"),
          path("locus_similarity_db.mmi"), emit: catalog
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def dedup_tolerance = params.quant_locus_dedup_tolerance ?: 10
    def chunk_size       = params.quant_chunk_size ?: 100
    def chunk_seed       = params.quant_chunk_seed ?: 42
    """
    quant_catalog_from_bed12.py --bed ${bed12} --tsv ${conf_tsv} --output ${unit_id}_catalog.tsv

    quant_locus_dedup.py \\
        --catalog ${unit_id}_catalog.tsv \\
        --tolerance ${dedup_tolerance} \\
        --prefix ${unit_id}

    quant_build_refs.py \\
        --catalog ${unit_id}_deduped_metadata.tsv \\
        --genome_fasta ${fasta} \\
        --outdir .

    quant_chunk_prep.py \\
        --catalog ${unit_id}_deduped_metadata.tsv \\
        --refs_fasta circle_refs.fa \\
        --refs_lengths circle_refs_lengths.tsv \\
        --chunk_size ${chunk_size} \\
        --seed ${chunk_seed} \\
        --outdir chunks

    quant_build_similarity_db.py \\
        --refs_fasta circle_refs.fa \\
        --refs_lengths circle_refs_lengths.tsv \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """
}
