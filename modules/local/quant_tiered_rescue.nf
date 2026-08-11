process QUANT_TIERED_RESCUE {
    tag "$meta.id"
    label 'process_high'
    container "${ (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/21497240/files/nanocirc-quant-latest.sif?download=1' :
        'quay.io/anrusakovich/nanocirc-quant@sha256:11534ee495f3717000a9b848600a395db092d22798a40eb4bccd533f7ca362ad' }"

    input:
    tuple val(meta), path(flagged_loci), path(tier1_counts), path(fastq),
          path(unit_refs_fasta), path(unit_refs_lengths), path(unit_similarity_db), path(unit_deduped_metadata)
    path fasta

    output:
    tuple val(meta), path("${meta.id}_quant_final_counts.tsv"), emit: counts
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quant_cross_locus_similarity.py \\
        --flagged_loci ${flagged_loci} \\
        --refs_fasta ${unit_refs_fasta} \\
        --refs_lengths ${unit_refs_lengths} \\
        --similarity_db ${unit_similarity_db} \\
        --deduped_metadata ${unit_deduped_metadata} \\
        --sample ${meta.id} \\
        --threads ${task.cpus}

    quant_tier2_rescue.py \\
        --flagged_similarity ${meta.id}_flagged_loci_similarity.tsv \\
        --deduped_metadata ${unit_deduped_metadata} \\
        --genome_fasta ${fasta} \\
        --reads ${fastq} \\
        --sample ${meta.id} \\
        --threads ${task.cpus}

    quant_tier3_rescue.py \\
        --flagged_similarity ${meta.id}_flagged_loci_similarity.tsv \\
        --deduped_metadata ${unit_deduped_metadata} \\
        --genome_fasta ${fasta} \\
        --reads ${fastq} \\
        --sample ${meta.id} \\
        --threads ${task.cpus}

    quant_finalize.py \\
        --deduped_metadata ${unit_deduped_metadata} \\
        --tier1_counts ${tier1_counts} \\
        --flagged_similarity ${meta.id}_flagged_loci_similarity.tsv \\
        --tier2_counts ${meta.id}_tier2_counts.tsv \\
        --tier3_counts ${meta.id}_tier3_counts.tsv \\
        --sample ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
