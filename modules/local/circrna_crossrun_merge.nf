process CIRCRNA_CROSSRUN_MERGE {
    tag "${meta.id}:${meta.tier}"
    label 'process_medium'

    container "${ (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container ?
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
    // cache_bust: smart_merge.py motif-informed pre-grouping + iso{N} label preservation on collision, 2026-08-24
    def sample_names = meta.sample_ids.join(' ')
    def n            = meta.sample_ids.size()
    def raw_prefix   = "${meta.id}_${meta.tier}_crossrun"
    prefix           = raw_prefix  // used in output block
    // A 1-sample group has nothing to corroborate against; requiring 2-sample
    // support (the normal balanced/high_confidence formula below) would drop
    // every locus. Floor both min_count and min_corroboration at 1 for n<2,
    // so a lone sample's own tier output passes through unchanged.
    def min_count    = n < 2 ? 1
                     : meta.tier == 'discovery' ? 1
                     : Math.max(2, Math.ceil(0.25 * n).toInteger())
    def min_corroboration = n < 2 ? 1 : params.crossrun_min_tool_agreement
    """
    # Only consensus_hybrid is computed here (see write_outputs()).
    # --conf_tsvs: cross-run mode, see cross_run_hybrid_entries().
    # --min_corroboration: min tool-agreement to keep a single-run structure.
    smart_merge.py \\
        --sample        ${raw_prefix} \\
        --tool_names    ${sample_names} \\
        --bed_files     ${bed_files} \\
        --conf_tsvs     ${conf_tsvs} \\
        --tolerance     ${params.circrna_bsj_tolerance} \\
        --struct_tolerance ${params.circrna_bsj_tolerance} \\
        --n_active      ${n} \\
        --min_corroboration ${min_corroboration} \\
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
    # PAIRS is empty for a 1-sample group (nothing to pair); omit --pairs then
    # rather than pass it an empty value.
    PAIRS_FLAG=()
    if [ "\${#PAIRS[@]}" -gt 0 ]; then
        PAIRS_FLAG=(--pairs "\${PAIRS[@]}")
    fi
    add_isoform_confidence.py \\
        --confidence  ${raw_prefix}_smart_consensus_hybrid_confidence.tsv \\
        "\${PAIRS_FLAG[@]}" \\
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
