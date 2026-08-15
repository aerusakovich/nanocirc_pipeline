/*
 * circRNA analysis subworkflow
 * Runs up to 4 detection tools, converts outputs to BED12,
 * and merges results when 2+ tools are active.
 */

include { ISOCIRC                } from '../../modules/local/isocirc'
include { CIRCFL_SEQ             } from '../../modules/local/circfl_seq'
include { PREPARE_GENOME         } from '../../modules/local/prepare_genome'
include { CIRI_LONG              } from '../../modules/local/ciri_long'
include { CIRCNICK_LRS           } from '../../modules/local/circnick_lrs'
include { CIRCNICK_LIFTOVER      } from '../../modules/local/circnick_liftover'
include { CIRILONG_TO_BED12      } from '../../modules/local/cirilong_to_bed12'
include { CIRCNICK_TO_BED12      } from '../../modules/local/circnick_to_bed12'
include { CIRCRNA_BEDTOOLS_PAIRS } from '../../modules/local/circrna_bedtools_pairs'
// Legacy benchmark modes, see misc/legacy_modules/
include { CIRCRNA_MERGE          } from '../../misc/legacy_modules/circrna_merge'
include { CIRCRNA_EXON_MERGE     } from '../../misc/legacy_modules/circrna_exon_merge'
include { CIRCRNA_SMART_MERGE    } from '../../modules/local/circrna_smart_merge'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_BALANCED           } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_HIGH_CONFIDENCE    } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_CONSENSUS_NO_LOW   } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_CONSENSUS_TRUSTED  } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_XSTRUCT_NO_LOW     } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_XSTRUCT_TRUSTED    } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_XSTRUCT_HIGH_ONLY  } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_PRIORITY_NO_LOW    } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_PRIORITY_TRUSTED   } from '../../modules/local/circrna_confidence_filter'
include { CIRCRNA_CONFIDENCE_FILTER as CIRCRNA_FILTER_PRIORITY_HIGH_ONLY } from '../../modules/local/circrna_confidence_filter'
include { GTF_TO_FEATURE_BED     } from '../../modules/local/gtf_to_feature_bed'
include { CIRCRNA_FINALIZE       } from '../../modules/local/circrna_finalize'
include { CIRCRNA_ANNOTATE       } from './circrna_annotate'

workflow CIRCRNA_ANALYSIS {

    take:
    ch_fastq              // channel: [ val(meta), path(fastq) ]
    fasta                 // path: reference genome FASTA
    gtf                   // path: gene annotation GTF
    circrna_db            // path: circRNA database

    main:

    def ch_versions = channel.empty()

    def runIsocirc       = Utils.asBool(params.run_isocirc)
    def runCircfl        = Utils.asBool(params.run_circfl)
    def runCirilong      = Utils.asBool(params.run_cirilong)
    def runCircnick      = Utils.asBool(params.run_circnick)
    def runBenchmarkModes = Utils.asBool(params.run_benchmark_modes)
    def skipAnnotation    = Utils.asBool(params.skip_annotation)

    // Gene/exon BED files derived from GTF, used for type classification
    GTF_TO_FEATURE_BED(gtf)
    ch_versions = ch_versions.mix(GTF_TO_FEATURE_BED.out.versions)
    def ch_gene_bed = GTF_TO_FEATURE_BED.out.gene_bed
    def ch_exon_bed = GTF_TO_FEATURE_BED.out.exon_bed

    // isocirc
    def ch_isocirc_bed = channel.empty()
    if (runIsocirc) {
        ISOCIRC ( ch_fastq, fasta, gtf, circrna_db )
        ch_isocirc_bed = ISOCIRC.out.bed
        ch_versions    = ch_versions.mix(ISOCIRC.out.versions.first())
    }

    // circFL-seq
    def ch_circfl_bed = channel.empty()
    if (runCircfl) {
        CIRCFL_SEQ ( ch_fastq, fasta, gtf )
        ch_circfl_bed = CIRCFL_SEQ.out.bed
        ch_versions   = ch_versions.mix(CIRCFL_SEQ.out.versions.first())
    }

    // CIRI-long
    def ch_cirilong_bed = channel.empty()
    if (runCirilong) {
        PREPARE_GENOME ( fasta )
        CIRI_LONG (
            ch_fastq,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.index,
            gtf,
            circrna_db
        )
        def ch_cirilong_info_isoforms = CIRI_LONG.out.info
            .map { m, f -> [m.id, m, f] }
            .join(CIRI_LONG.out.isoforms.map { m, f -> [m.id, f] })
            .map { _id, m, info, isoforms -> [m, info, isoforms] }
        CIRILONG_TO_BED12 ( ch_cirilong_info_isoforms )
        ch_cirilong_bed = CIRILONG_TO_BED12.out.bed12
        ch_versions     = ch_versions.mix(PREPARE_GENOME.out.versions)
        ch_versions     = ch_versions.mix(CIRI_LONG.out.versions.first())
        ch_versions     = ch_versions.mix(CIRILONG_TO_BED12.out.versions.first())
    }

    // circnick-lrs
    def ch_circnick_bed = channel.empty()
    if (runCircnick) {
        CIRCNICK_LRS ( ch_fastq, params.circnick_species )
        ch_versions = ch_versions.mix(CIRCNICK_LRS.out.versions.first())

        if (params.circnick_liftover_chain) {
            def ch_chain = channel.fromPath(params.circnick_liftover_chain, checkIfExists: true)

            CIRCNICK_LIFTOVER (
                CIRCNICK_LRS.out.annotated,
                CIRCNICK_LRS.out.exon_usage,
                CIRCNICK_LRS.out.intron_cov,
                ch_chain.first()
            )
            ch_versions = ch_versions.mix(CIRCNICK_LIFTOVER.out.versions.first())

            CIRCNICK_TO_BED12 (
                CIRCNICK_LIFTOVER.out.annotated,
                CIRCNICK_LIFTOVER.out.exon_usage,
                CIRCNICK_LIFTOVER.out.intron_cov
            )
        } else {
            CIRCNICK_TO_BED12 (
                CIRCNICK_LRS.out.annotated,
                CIRCNICK_LRS.out.exon_usage,
                CIRCNICK_LRS.out.intron_cov
            )
        }

        ch_circnick_bed = CIRCNICK_TO_BED12.out.bed12
        ch_versions     = ch_versions.mix(CIRCNICK_TO_BED12.out.versions.first())
    }

    // Native per-tool outputs (not BED12) needed by the quantification subworkflow's
    // tool_readset_union step.
    def ch_isocirc_native  = runIsocirc  ? ISOCIRC.out.expr
                                         : ch_fastq.map { m, _f -> [m, file('NO_FILE_ISOCIRC')] }
    def ch_circfl_native   = runCircfl   ? CIRCFL_SEQ.out.expr.map { m, files ->
                                  def fl = files instanceof List ? files : [files]
                                  [m, fl.find { f -> f.toString().contains('/mRG/') } ?: fl[0]]
                              }          : ch_fastq.map { m, _f -> [m, file('NO_FILE_CIRCFL')] }

    def ch_cirilong_reads  = runCirilong
        ? ch_fastq.map { m, _f -> [m.id, m] }
            .join(CIRI_LONG.out.reads.map { m, f -> [m.id, f] }, remainder: true)
            .map { _id, m, f -> [m, f ?: file('NO_FILE_CIRILONG')] }
        : ch_fastq.map { m, _f -> [m, file('NO_FILE_CIRILONG')] }
    def ch_circnick_native = runCircnick ? CIRCNICK_LRS.out.annotated
                                         : ch_fastq.map { m, _f -> [m, file('NO_FILE_CIRCNICK')] }

    def ch_all_beds = channel.empty()

    if (runIsocirc) {
        ch_all_beds = ch_all_beds.mix(
            ch_isocirc_bed.map { meta, bed -> [ meta, 'isocirc', bed ] }
        )
    }
    if (runCircfl) {
        ch_all_beds = ch_all_beds.mix(
            ch_circfl_bed.map { meta, bed -> [ meta, 'circfl', bed ] }
        )
    }
    if (runCirilong) {
        ch_all_beds = ch_all_beds.mix(
            ch_cirilong_bed.map { meta, bed -> [ meta, 'cirilong', bed ] }
        )
    }
    if (runCircnick) {
        ch_all_beds = ch_all_beds.mix(
            ch_circnick_bed.map { meta, bed -> [ meta, 'circnick', bed ] }
        )
    }

    // Group by sample: [ meta, [tool_names], [bed_files] ]
    def ch_beds_collected = ch_all_beds
        .groupTuple()
        .map { meta, tool_names, bed_files ->
            [ meta, tool_names, bed_files ]
        }

    def ch_n_active = runIsocirc  ? 1 : 0
    ch_n_active = runCircfl   ? ch_n_active + 1 : ch_n_active
    ch_n_active = runCirilong ? ch_n_active + 1 : ch_n_active
    ch_n_active = runCircnick ? ch_n_active + 1 : ch_n_active

    if (ch_n_active < 4) {
        log.warn(
            "[${workflow.manifest.name}] Only ${ch_n_active} of 4 detection tools are active. " +
            "Tool consensus scores (Low/Medium/High) reflect agreement among the " +
            "tools that ran. A 'High' score from ${ch_n_active} tools is not " +
            "equivalent to 'High' from all 4 tools."
        )
    }

    if (ch_n_active >= 1) {

        CIRCRNA_BEDTOOLS_PAIRS ( ch_beds_collected )
        ch_versions = ch_versions.mix(CIRCRNA_BEDTOOLS_PAIRS.out.versions.first())

        def ch_for_merge = ch_beds_collected
            .join(CIRCRNA_BEDTOOLS_PAIRS.out.pairs, by: 0, remainder: true)
            .map { meta, tool_names, bed_files, pairs -> [ meta, tool_names, bed_files, pairs ?: [] ] }

        // Legacy benchmark modes, see misc/legacy_modules/
        if (runBenchmarkModes) {
            CIRCRNA_MERGE (
                ch_for_merge.map { meta, tool_names, bed_files, pairs -> [ meta, tool_names, bed_files ] },
                ch_for_merge.map { meta, tool_names, bed_files, pairs -> [ meta, pairs ] },
                ch_n_active
            )
            ch_versions = ch_versions.mix(CIRCRNA_MERGE.out.versions.first())

            CIRCRNA_EXON_MERGE (
                CIRCRNA_BEDTOOLS_PAIRS.out.pairs,
                ch_n_active
            )
            ch_versions = ch_versions.mix(CIRCRNA_EXON_MERGE.out.versions.first())
        }

        CIRCRNA_SMART_MERGE (
            ch_for_merge.map { meta, tool_names, bed_files, pairs -> [ meta, tool_names, bed_files ] },
            ch_for_merge.map { meta, tool_names, bed_files, pairs -> [ meta, pairs ] },
            ch_n_active
        )
        ch_versions = ch_versions.mix(CIRCRNA_SMART_MERGE.out.versions.first())

        def ch_hybrid = CIRCRNA_SMART_MERGE.out.hybrid_bed
            .join(CIRCRNA_SMART_MERGE.out.hybrid_conf, by: 0)

        def ch_for_balanced        = ch_hybrid.map { meta, bed, tsv -> [ meta + [category: 'hybrid'], bed, tsv ] }
        def ch_for_high_confidence = ch_hybrid.map { meta, bed, tsv -> [ meta + [category: 'hybrid'], bed, tsv ] }
        def ch_for_consensus_trusted = ch_hybrid.map { meta, bed, tsv -> [ meta + [category: 'hybrid'], bed, tsv ] }

        CIRCRNA_FILTER_BALANCED          ( ch_for_balanced )
        CIRCRNA_FILTER_HIGH_CONFIDENCE   ( ch_for_high_confidence )
        CIRCRNA_FILTER_CONSENSUS_TRUSTED ( ch_for_consensus_trusted.map { meta, bed, tsv -> [ meta + [category: 'balanced_recall'], bed, tsv ] } )
        ch_versions = ch_versions.mix(CIRCRNA_FILTER_BALANCED.out.versions.first())
        ch_versions = ch_versions.mix(CIRCRNA_FILTER_CONSENSUS_TRUSTED.out.versions.first())

        // Only used below by the benchmark-mode smart_consensus_no_low diagnostic.
        def ch_consensus = CIRCRNA_SMART_MERGE.out.consensus_bed
            .join(CIRCRNA_SMART_MERGE.out.consensus_conf, by: 0)

        if (runBenchmarkModes) {
            def ch_xstruct = CIRCRNA_SMART_MERGE.out.consensus_xstruct_bed
                .join(CIRCRNA_SMART_MERGE.out.consensus_xstruct_conf, by: 0)

            def ch_for_priority = CIRCRNA_SMART_MERGE.out.priority_bed
                .join(CIRCRNA_SMART_MERGE.out.priority_conf, by: 0)

            CIRCRNA_FILTER_CONSENSUS_NO_LOW  ( ch_consensus.map { meta, bed, tsv -> [ meta + [category: 'smart_consensus_no_low'],      bed, tsv ] } )
            CIRCRNA_FILTER_XSTRUCT_NO_LOW    ( ch_xstruct.map   { meta, bed, tsv -> [ meta + [category: 'smart_consensus_xstruct_no_low'],  bed, tsv ] } )
            CIRCRNA_FILTER_XSTRUCT_TRUSTED   ( ch_xstruct.map   { meta, bed, tsv -> [ meta + [category: 'smart_consensus_xstruct_filtered'], bed, tsv ] } )
            CIRCRNA_FILTER_XSTRUCT_HIGH_ONLY ( ch_xstruct.map   { meta, bed, tsv -> [ meta + [category: 'smart_consensus_xstruct_high_only'], bed, tsv ] } )
            CIRCRNA_FILTER_PRIORITY_NO_LOW   ( ch_for_priority.map { meta, bed, tsv -> [ meta + [category: 'smart_priority_no_low'],    bed, tsv ] } )
            CIRCRNA_FILTER_PRIORITY_TRUSTED  ( ch_for_priority.map { meta, bed, tsv -> [ meta + [category: 'smart_priority_filtered'],  bed, tsv ] } )
            CIRCRNA_FILTER_PRIORITY_HIGH_ONLY ( ch_for_priority.map { meta, bed, tsv -> [ meta + [category: 'smart_priority_high_only'], bed, tsv ] } )
        }

        // Discovery is hybrid emitted directly (no filter applied)
        def ch_discovery_for_annotate = ch_hybrid
            .map { meta, bed, tsv -> [ meta + [category: 'discovery'], bed, tsv ] }

        if (!skipAnnotation) {

            def ch_for_annotate = ch_discovery_for_annotate
                .mix( CIRCRNA_FILTER_BALANCED.out.bed          .join(CIRCRNA_FILTER_BALANCED.out.conf,          by: 0).map { meta, bed, tsv -> [ meta + [category: 'balanced_precision'], bed, tsv ] } )
                .mix( CIRCRNA_FILTER_CONSENSUS_TRUSTED.out.bed .join(CIRCRNA_FILTER_CONSENSUS_TRUSTED.out.conf, by: 0).map { meta, bed, tsv -> [ meta + [category: 'balanced_recall'],    bed, tsv ] } )
                .mix( CIRCRNA_FILTER_HIGH_CONFIDENCE.out.bed   .join(CIRCRNA_FILTER_HIGH_CONFIDENCE.out.conf,   by: 0).map { meta, bed, tsv -> [ meta + [category: 'high_confidence'],    bed, tsv ] } )

            if (runBenchmarkModes) {
                ch_for_annotate = ch_for_annotate
                    .mix( CIRCRNA_MERGE.out.strict_union_bed  .join(CIRCRNA_MERGE.out.strict_union_conf,  by: 0).map { meta, bed, tsv -> [ meta + [category: 'strict_union'],  bed, tsv ] } )
                    .mix( CIRCRNA_MERGE.out.strict_inter_bed  .join(CIRCRNA_MERGE.out.strict_inter_conf,  by: 0).map { meta, bed, tsv -> [ meta + [category: 'strict_inter'],  bed, tsv ] } )
                    .mix( CIRCRNA_MERGE.out.relaxed_union_bed .join(CIRCRNA_MERGE.out.relaxed_union_conf, by: 0).map { meta, bed, tsv -> [ meta + [category: 'relaxed_union'], bed, tsv ] } )
                    .mix( CIRCRNA_MERGE.out.relaxed_inter_bed .join(CIRCRNA_MERGE.out.relaxed_inter_conf, by: 0).map { meta, bed, tsv -> [ meta + [category: 'relaxed_inter'], bed, tsv ] } )
                    .mix( CIRCRNA_EXON_MERGE.out.exon_union_bed.join(CIRCRNA_EXON_MERGE.out.exon_union_conf, by: 0).map { meta, bed, tsv -> [ meta + [category: 'exon_union'], bed, tsv ] } )
                    .mix( CIRCRNA_EXON_MERGE.out.exon_inter_bed.join(CIRCRNA_EXON_MERGE.out.exon_inter_conf, by: 0).map { meta, bed, tsv -> [ meta + [category: 'exon_inter'], bed, tsv ] } )
            }

            CIRCRNA_ANNOTATE(ch_for_annotate, fasta, gtf)
            // CIRCRNA_ANNOTATE uses nf-core modules with `topic: versions`.
            // Versions are collected automatically, no ch_versions mixing needed.

            // ── Finalize: type + expression + clean TSV ─────────────────────────────
            // annotated_tsv: N items/sample. expr channels: 1/sample. combine(by:0), not join().
            // circfl: pick mRG pass file when present.
            // Placeholder file names must be distinct per tool: Nextflow errors on
            // duplicate staged filenames when 2+ tools are inactive at once.
            def ch_iso_expr  = runIsocirc  ? ISOCIRC.out.expr.map { m, f -> [m.id, f] }
                                           : ch_fastq.map { m, _f -> [m.id, file('NO_FILE_ISOCIRC')] }
            def ch_fl_expr   = runCircfl   ? CIRCFL_SEQ.out.expr.map { m, files ->
                                   def fl = files instanceof List ? files : [files]
                                   [m.id, fl.find { f -> f.toString().contains('/mRG/') } ?: fl[0]]
                               }           : ch_fastq.map { m, _f -> [m.id, file('NO_FILE_CIRCFL')] }
            def ch_nick_expr = runCircnick ? CIRCNICK_LRS.out.annotated.map { m, f -> [m.id, f] }
                                           : ch_fastq.map { m, _f -> [m.id, file('NO_FILE_CIRCNICK')] }
            // ciri_long: guarantee one entry per sample even when expression file is absent
            def ch_ciri_expr = runCirilong
                ? ch_fastq.map { m, _f -> [m.id] }
                    .join(CIRI_LONG.out.expr.map { m, f -> [m.id, f] }, remainder: true)
                    .map { id, f -> [id, f ?: file('NO_FILE_CIRILONG')] }
                : ch_fastq.map { m, _f -> [m.id, file('NO_FILE_CIRILONG')] }

            def ch_for_finalize = CIRCRNA_ANNOTATE.out.annotated_tsv
                .map     { meta, tsv -> [meta.id, meta, tsv] }
                .combine ( ch_iso_expr,  by: 0 )
                .combine ( ch_fl_expr,   by: 0 )
                .combine ( ch_nick_expr, by: 0 )
                .combine ( ch_ciri_expr, by: 0 )
                .map     { id, meta, tsv, iso, fl, nick, ciri -> [meta, tsv, iso, ciri, fl, nick] }

            CIRCRNA_FINALIZE(ch_for_finalize, ch_gene_bed, ch_exon_bed)
            ch_versions = ch_versions.mix(CIRCRNA_FINALIZE.out.versions.first())
        }
    }

    emit:
    // Per-tool BED12 outputs (always present when tool is active)
    isocirc_bed12   = ch_isocirc_bed
    circfl_bed12    = ch_circfl_bed
    cirilong_bed12  = ch_cirilong_bed
    circnick_bed12  = ch_circnick_bed

    // Feature BED files (for type classification, passed to crossrun merge)
    gene_bed = ch_gene_bed
    exon_bed = ch_exon_bed

    // Native per-tool outputs (for the quantification subworkflow's tool_readset_union)
    isocirc_native  = ch_isocirc_native
    circfl_native   = ch_circfl_native
    cirilong_reads  = ch_cirilong_reads
    circnick_native = ch_circnick_native

    // Discovery: hybrid, all isoforms, unfiltered. With 1 tool active, every
    // score is trivially "100% of active tools agree" (see add_isoform_confidence.py's
    // count_to_score), so all 4 tiers below end up identical to this: there is no
    // second tool to disagree with, so nothing is ever Low/Medium.
    discovery_bed  = ch_n_active >= 1 ? CIRCRNA_SMART_MERGE.out.hybrid_bed  : channel.empty()
    discovery_conf = ch_n_active >= 1 ? CIRCRNA_SMART_MERGE.out.hybrid_conf : channel.empty()
    // Balanced (precision): hybrid + no_low filter
    balanced_bed   = ch_n_active >= 1 ? CIRCRNA_FILTER_BALANCED.out.bed        : channel.empty()
    balanced_conf  = ch_n_active >= 1 ? CIRCRNA_FILTER_BALANCED.out.conf       : channel.empty()
    // Balanced (recall): hybrid + trusted_only filter
    balanced_recall_bed  = ch_n_active >= 1 ? CIRCRNA_FILTER_CONSENSUS_TRUSTED.out.bed  : channel.empty()
    balanced_recall_conf = ch_n_active >= 1 ? CIRCRNA_FILTER_CONSENSUS_TRUSTED.out.conf : channel.empty()
    // High-confidence: hybrid + high_only filter
    high_conf_bed  = ch_n_active >= 1 ? CIRCRNA_FILTER_HIGH_CONFIDENCE.out.bed  : channel.empty()
    high_conf_conf = ch_n_active >= 1 ? CIRCRNA_FILTER_HIGH_CONFIDENCE.out.conf : channel.empty()

    // Legacy merge modes (--run_benchmark_modes) are strict/relaxed/exon
    // union+intersection across >=2 tools' raw calls; not meaningful with 1
    // tool, so these stay off for a single-tool run even with the flag set.
    strict_union_bed    = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_MERGE.out.strict_union_bed    : channel.empty()
    strict_union_conf   = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_MERGE.out.strict_union_conf   : channel.empty()
    strict_inter_bed    = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_MERGE.out.strict_inter_bed    : channel.empty()
    strict_inter_conf   = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_MERGE.out.strict_inter_conf   : channel.empty()
    relaxed_union_bed   = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_MERGE.out.relaxed_union_bed   : channel.empty()
    relaxed_union_conf  = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_MERGE.out.relaxed_union_conf  : channel.empty()
    relaxed_inter_bed   = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_MERGE.out.relaxed_inter_bed   : channel.empty()
    relaxed_inter_conf  = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_MERGE.out.relaxed_inter_conf  : channel.empty()
    exon_union_bed      = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_EXON_MERGE.out.exon_union_bed  : channel.empty()
    exon_union_conf     = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_EXON_MERGE.out.exon_union_conf : channel.empty()
    exon_inter_bed      = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_EXON_MERGE.out.exon_inter_bed  : channel.empty()
    exon_inter_conf     = (ch_n_active >= 2 && runBenchmarkModes) ? CIRCRNA_EXON_MERGE.out.exon_inter_conf : channel.empty()

    annotated_gff = (!skipAnnotation && ch_n_active >= 1) ? CIRCRNA_ANNOTATE.out.annotated_gff : channel.empty()
    spliced_fasta = (!skipAnnotation && ch_n_active >= 1) ? CIRCRNA_ANNOTATE.out.spliced_fasta : channel.empty()
    annotated_tsv = (!skipAnnotation && ch_n_active >= 1) ? CIRCRNA_ANNOTATE.out.annotated_tsv : channel.empty()
    clean_tsv     = (!skipAnnotation && ch_n_active >= 1) ? CIRCRNA_FINALIZE.out.clean          : channel.empty()

    versions = ch_versions
}
