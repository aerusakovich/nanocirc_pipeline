/*
 * circRNA quantification subworkflow
 *
 * Remap-based circRNA quantification: tier1 (chunked remap and classify
 * against synthetic circle references), then tier1.5 (overlap-cluster
 * rescue), tier2 (targeted low-coverage rescue), tier3 (gene-family
 * rescue), and finalize. Ground truth, structural matching and plots are
 * simulator-only benchmarking code and are not part of this subworkflow.
 *
 * Catalog sharing: circle_refs.fa, chunks and locus_similarity_db are
 * built once per dataset, not once per run. Only the per-run counts are
 * separate.
 *   - Per-sample quantification (no crossrun grouping): each sample is
 *     quantified against its own discovery catalog.
 *   - Per-group quantification (--run_crossrun_merge): every run in a
 *     group shares one catalog built from that group's crossrun
 *     discovery output. This is needed for a coherent per-locus count
 *     matrix across samples.
 * The caller (workflows/nanocirc.nf) decides which mode applies and
 * builds `ch_unit_catalog` (one catalog per "quant unit", a group or a
 * singleton sample) and `ch_run_to_unit` (maps every actual sample/run to
 * its unit). This subworkflow does not care which mode produced them.
 *
 * Only the reference-genome minimap2 index (quant_genome_index, no
 * catalog dependency at all) is always shared across every unit and run.
 */

include { QUANT_GENOME_INDEX    } from '../../modules/local/quant_genome_index'
include { QUANT_PREPARE_CATALOG } from '../../modules/local/quant_prepare_catalog'
include { QUANT_TIER1           } from '../../modules/local/quant_tier1'
include { QUANT_COVERAGE_FLAG   } from '../../modules/local/quant_coverage_flag'
include { QUANT_TIERED_RESCUE   } from '../../modules/local/quant_tiered_rescue'

workflow CIRCRNA_QUANTIFY {

    take:
    ch_fastq             // channel: [ val(meta), path(fastq) ], full raw reads per sample
    ch_sample_discovery  // channel: [ val(meta), path(bed12), path(conf_tsv) ], each sample's own discovery catalog (always per-sample)
    ch_unit_catalog      // channel: [ val(unit_id), path(bed12), path(conf_tsv) ], one per quant unit (group or singleton sample)
    ch_run_to_unit       // channel: [ val(meta), val(unit_id) ], maps every sample/run to its quant unit
    fasta                // path: reference genome FASTA
    ch_isocirc_native    // channel: [ val(meta), path ], native isocirc.out (or NO_FILE_ISOCIRC)
    ch_circfl_native     // channel: [ val(meta), path ], native circFL_Normal_pass.txt (or NO_FILE_CIRCFL)
    ch_cirilong_reads    // channel: [ val(meta), path ], native CIRI-long .reads (or NO_FILE_CIRILONG)
    ch_circnick_native   // channel: [ val(meta), path ], native circnick annotated.txt (or NO_FILE_CIRCNICK)

    main:

    def ch_versions = channel.empty()

    // Reference-genome index, shared and broadcast to every sample.
    QUANT_GENOME_INDEX(fasta)
    ch_versions = ch_versions.mix(QUANT_GENOME_INDEX.out.versions)
    def ch_genome_index = QUANT_GENOME_INDEX.out.index.first()

    // ── Catalog build (once per quant unit: group or singleton sample) ──────────
    QUANT_PREPARE_CATALOG(ch_unit_catalog, fasta)
    ch_versions = ch_versions.mix(QUANT_PREPARE_CATALOG.out.versions.first())

    // Broadcast each unit's catalog artifacts to every run that belongs to it
    def ch_run_catalog = ch_run_to_unit
        .map { meta, unit_id -> [unit_id, meta] }
        .combine(QUANT_PREPARE_CATALOG.out.catalog, by: 0)
        .map { _unit_id, meta, deduped_metadata, cluster_map, refs_fasta, refs_lengths, chunks_dir, similarity_db ->
            [meta, deduped_metadata, cluster_map, refs_fasta, refs_lengths, chunks_dir, similarity_db]
        }

    // ── Tier1 + Tier1.5: chunked remap-and-classify + overlap-cluster rescue ────
    QUANT_TIER1(
        ch_fastq
            .join(ch_run_catalog, by: 0)
            .map { meta, fastq, dedup, _cmap, rfasta, rlen, chunks, _simdb ->
                [meta, fastq, dedup, rfasta, rlen, chunks]
            },
        ch_genome_index
    )
    ch_versions = ch_versions.mix(QUANT_TIER1.out.versions.first())

    // ── Old-tool-count baseline + coverage flagging ──────────────────────────────
    QUANT_COVERAGE_FLAG(
        ch_sample_discovery
            .join(ch_isocirc_native, by: 0)
            .join(ch_circfl_native, by: 0)
            .join(ch_cirilong_reads, by: 0)
            .join(ch_circnick_native, by: 0)
            .join(ch_run_catalog.map { meta, dedup, cmap, _rfasta, _rlen, _chunks, _simdb -> [meta, dedup, cmap] }, by: 0)
            .join(QUANT_TIER1.out.counts, by: 0)
    )
    ch_versions = ch_versions.mix(QUANT_COVERAGE_FLAG.out.versions.first())

    // ── Gene-family similarity check + Tier2 (targeted) + Tier3 (gene-family) ───
    // rescue + finalize ──────────────────────────────────────────────────────────
    QUANT_TIERED_RESCUE(
        QUANT_COVERAGE_FLAG.out.flagged
            .join(QUANT_TIER1.out.counts, by: 0)
            .join(ch_fastq, by: 0)
            .join(ch_run_catalog.map { meta, dedup, _cmap, rfasta, rlen, _chunks, simdb -> [meta, rfasta, rlen, simdb, dedup] }, by: 0),
        fasta
    )
    ch_versions = ch_versions.mix(QUANT_TIERED_RESCUE.out.versions.first())

    emit:
    final_counts = QUANT_TIERED_RESCUE.out.counts   // channel: [ val(meta), path(quant_final_counts.tsv) ]
    // Must join onto clean.tsv via cluster_map, not directly: skipping this
    // hop silently zeroes ~19% of loci's counts (see quant_append_counts.py).
    cluster_map  = ch_run_catalog.map { meta, _dedup, cmap, _rfasta, _rlen, _chunks, _simdb -> [meta, cmap] }
    versions     = ch_versions
}
