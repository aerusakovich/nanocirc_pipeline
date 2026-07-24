process QUANT_TIER1 {
    tag "$meta.id"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://zenodo.org/records/21497240/files/nanocirc-quant-latest.sif?download=1' :
        'quay.io/anrusakovich/nanocirc-quant:latest' }"

    input:
    tuple val(meta), path(fastq), path(deduped_metadata), path(refs_fasta), path(refs_lengths), path(chunks_dir)
    path genome_index

    output:
    tuple val(meta), path("${meta.id}_tier1p5_counts.tsv"), emit: counts
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quant_genome_align.py \\
        --genome_index ${genome_index} \\
        --reads ${fastq} \\
        --sample ${meta.id} \\
        --threads ${task.cpus}

    quant_chunk_remap.py \\
        --chunks_dir ${chunks_dir} \\
        --genome_bam ${meta.id}.genome_align.sorted.bam \\
        --sample ${meta.id} \\
        --threads ${task.cpus}

    quant_overlap_rescue.py \\
        --catalog ${deduped_metadata} \\
        --tier1_counts ${meta.id}_tier1_counts.tsv \\
        --genome_bam ${meta.id}.genome_align.sorted.bam \\
        --refs_fasta ${refs_fasta} \\
        --refs_lengths ${refs_lengths} \\
        --sample ${meta.id} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
