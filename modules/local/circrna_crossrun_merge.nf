process CIRCRNA_CROSSRUN_MERGE {
    tag "${meta.id}:${meta.tier}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybedtools:0.12.0--py39h475c85d_0' :
        'quay.io/biocontainers/pybedtools:0.12.0--py39h475c85d_0' }"

    input:
    tuple val(meta), path(bed_files), path(conf_tsvs)
    path gene_bed
    path exon_bed

    output:
    tuple val(meta), path("${prefix}.bed12"),            emit: bed
    tuple val(meta), path("${prefix}_confidence.tsv"),   emit: confidence
    tuple val(meta), path("${prefix}_clean.tsv"),        emit: clean
    path  "versions.yml",                                emit: versions

    script:
    def sample_names = meta.sample_ids.join(' ')
    def n            = meta.sample_ids.size()
    def raw_prefix   = "${meta.id}_${meta.tier}_crossrun"
    prefix           = raw_prefix  // used in output block
    def min_count    = meta.tier == 'discovery'                                    ? 1
                     : (meta.tier == 'balanced_precision' || meta.tier == 'balanced_recall') ? Math.max(2, Math.ceil(0.25 * n).toInteger())
                     : /* high_confidence */                                        Math.ceil(0.75 * n).toInteger()
    """
    # smart_merge: makes all 4 mode files, only the hybrid ones are used below.
    # The consensus, xstruct and priority files stay in the work dir, unused.
    #
    # --conf_tsvs turns on cross-run merge mode. Here each sample is
    # already a merged catalog, not one raw call per tool, so its bed12 can
    # have several records (main plus isoN) at the same locus. Without
    # --conf_tsvs, smart_merge.py reduces each sample to one record by raw
    # score before voting. That can drop a sample's own correct 'main' pick
    # if another record from the same sample has a higher score.
    # --conf_tsvs instead groups every record by its real structure and
    # ranks groups by summed struct_agree_count (total tool agreement).
    # See cross_run_hybrid_entries() in smart_merge.py.
    #
    # --min_corroboration drops a structure seen in only 1 run, unless that
    # run's own tool agreement meets this number. Tests showed this cuts
    # cross-run false positives by about 34%, for a small recall cost.
    smart_merge.py \\
        --sample        ${raw_prefix} \\
        --tool_names    ${sample_names} \\
        --bed_files     ${bed_files} \\
        --conf_tsvs     ${conf_tsvs} \\
        --tolerance     ${params.circrna_bsj_tolerance} \\
        --struct_tolerance ${params.circrna_bsj_tolerance} \\
        --n_active      ${n} \\
        --min_corroboration ${params.crossrun_min_tool_agreement} \\
        --outdir        .

    # ── Pairwise bedtools intersect ────────────────────────────────────────────
    BEDS=(${bed_files})
    NAMES=(${sample_names})
    PAIRS=()
    for (( i=0; i<${n}; i++ )); do
        for (( j=i+1; j<${n}; j++ )); do
            pf="\${NAMES[i]}_vs_\${NAMES[j]}.pairs.bed"
            bedtools intersect \\
                -a "\${BEDS[i]}" -b "\${BEDS[j]}" -split -wo > "\${pf}"
            PAIRS+=("\${pf}")
        done
    done

    # ── Add isoform confidence ─────────────────────────────────────────────────
    add_isoform_confidence.py \\
        --confidence  ${raw_prefix}_smart_consensus_hybrid_confidence.tsv \\
        --pairs       "\${PAIRS[@]}" \\
        --min_overlap ${params.circrna_isoform_overlap} \\
        --n_active    ${n} \\
        --strip_isoform_suffix \\
        --output      ${raw_prefix}_smart_consensus_hybrid_confidence.tsv

    # ── Annotate and produce final outputs ────────────────────────────────────
    # Reads the raw hybrid TSV (above) without modifying it;
    # writes prefix_confidence.tsv (full intermediate), prefix_clean.tsv and prefix.bed12.
    crossrun_annotate.py \\
        --input_tsv    ${raw_prefix}_smart_consensus_hybrid_confidence.tsv \\
        --input_bed    ${raw_prefix}_smart_consensus_hybrid.bed12 \\
        --sample_names ${sample_names} \\
        --sample_tsvs  ${conf_tsvs} \\
        --gene_bed     ${gene_bed} \\
        --exon_bed     ${exon_bed} \\
        --bsj_tol      ${params.circrna_bsj_tolerance} \\
        --min_count    ${min_count} \\
        --prefix       ${raw_prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | head -1 | sed 's/bedtools v//')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
