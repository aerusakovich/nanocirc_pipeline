/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { NANOPLOT               } from '../modules/nf-core/nanoplot/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { FASTQ_RENAME           } from '../modules/local/fastq_rename'
include { CIRCRNA_ANALYSIS            } from '../subworkflows/local/circrna_analysis'
include { CIRCRNA_QUANTIFY        } from '../subworkflows/local/circrna_quantify'
include { QUANT_APPEND_COUNTS     } from '../modules/local/quant_append_counts'
include { CIRCRNA_CROSSRUN_MERGE  } from '../modules/local/circrna_crossrun_merge'
include { FILTER_CONFIDENT_DISCOVERY } from '../modules/local/filter_confident_discovery'
include { BUILD_DESEQ2_MATRIX     } from '../modules/local/build_deseq2_matrix'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanocirc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOCIRC {

    take:
    ch_samplesheet // channel: [ val(meta), path(fastq) ] from --input samplesheet

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Coerce every boolean CLI flag used below once, up front. See Utils.asBool().
    def skipQc          = Utils.asBool(params.skip_qc)
    def skipFastqc      = Utils.asBool(params.skip_fastqc)
    def skipNanoplot    = Utils.asBool(params.skip_nanoplot)
    def skipMultiqc     = Utils.asBool(params.skip_multiqc)
    def runIsocirc      = Utils.asBool(params.run_isocirc)
    def runCircfl       = Utils.asBool(params.run_circfl)
    def runCirilong     = Utils.asBool(params.run_cirilong)
    def runCircnick     = Utils.asBool(params.run_circnick)
    def runCrossrunMerge = Utils.asBool(params.run_crossrun_merge)
    def runQuantify      = Utils.asBool(params.run_quantify)

    // Ensure all FASTQ files are named *.fastq.gz. NanoPlot needs this,
    // it detects input type from the extension.
    FASTQ_RENAME ( ch_samplesheet )
    ch_fastq = FASTQ_RENAME.out.fastq

    //
    // QC: FastQC and NanoPlot
    //
    if (!skipQc) {
        if (!skipFastqc) {
            FASTQC ( ch_fastq )
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ v -> v[1] })
            ch_versions      = ch_versions.mix(FASTQC.out.versions.first())
        }
        if (!skipNanoplot) {
            NANOPLOT ( ch_fastq )
            ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())
        }
    }

    //
    // SUBWORKFLOW: circRNA detection, BED12 conversion, merge and confidence scoring
    //
    if (!params.fasta) {
        error("CircRNA analysis requires '--fasta' (reference genome).")
    }
    if (!params.gtf) {
        error("CircRNA analysis requires '--gtf' (gene annotation).")
    }
    if ((runIsocirc || runCirilong) && !params.circrna_db) {
        error("isocirc and ciri-long require '--circrna_db' (circRNA database).")
    }
    if (!runIsocirc && !runCircfl && !runCirilong && !runCircnick) {
        error("At least one tool must be active. Use --run_isocirc, --run_circfl, --run_cirilong, or --run_circnick.")
    }
    if (runCircnick) {
        if (!params.circnick_species) {
            error("Parameter '--circnick_species' is required when '--run_circnick' is set. Valid options: 'mouse', 'human'")
        }
        if (params.circnick_species != 'mouse' && params.circnick_species != 'human') {
            error("Invalid --circnick_species: '${params.circnick_species}'. Valid options: 'mouse', 'human'")
        }
        if (!params.circnick_liftover_chain) {
            log.warn "No --circnick_liftover_chain provided. circnick-lrs will use built-in ${params.circnick_species == 'mouse' ? 'mm10' : 'hg19'} coordinates."
        }
    }

    if (runCrossrunMerge) {
        // Validate: every sample must declare a group when cross-run merge is on
        ch_samplesheet
            .filter { meta, _fastq -> !meta.group }
            .subscribe { meta, _fastq ->
                error("--run_crossrun_merge requires a 'group' column in the samplesheet, but sample '${meta.id}' has none.")
            }
    }

    CIRCRNA_ANALYSIS (
        ch_fastq,
        file(params.fasta, checkIfExists: true),
        file(params.gtf,   checkIfExists: true),
        params.circrna_db ? file(params.circrna_db, checkIfExists: true) : file('NO_FILE')
    )
    ch_versions = ch_versions.mix(CIRCRNA_ANALYSIS.out.versions)

    //
    // Cross-sample merge: group per-sample outputs by samplesheet 'group' field
    //
    if (runCrossrunMerge) {
        def group_by_tier = { bed_ch, conf_ch, tier_name ->
            bed_ch
                .filter  { meta, _bed -> meta.group }
                .join    ( conf_ch, by: 0 )
                .map     { meta, bed, tsv -> tuple(meta.group, meta.id, bed, tsv) }
                .groupTuple(sort: { entry -> entry[0] })
                .map     { grp, sample_ids, beds, tsvs ->
                    tuple([id: grp, sample_ids: sample_ids, tier: tier_name], beds, tsvs)
                }
        }

        def ch_crossrun = group_by_tier.call(
                CIRCRNA_ANALYSIS.out.discovery_bed,  CIRCRNA_ANALYSIS.out.discovery_conf,  'discovery')
            .mix(group_by_tier.call(
                CIRCRNA_ANALYSIS.out.balanced_bed,   CIRCRNA_ANALYSIS.out.balanced_conf,   'balanced_precision'))
            .mix(group_by_tier.call(
                CIRCRNA_ANALYSIS.out.balanced_recall_bed, CIRCRNA_ANALYSIS.out.balanced_recall_conf, 'balanced_recall'))
            .mix(group_by_tier.call(
                CIRCRNA_ANALYSIS.out.high_conf_bed,  CIRCRNA_ANALYSIS.out.high_conf_conf,  'high_confidence'))

        CIRCRNA_CROSSRUN_MERGE (
            ch_crossrun,
            CIRCRNA_ANALYSIS.out.gene_bed,
            CIRCRNA_ANALYSIS.out.exon_bed
        )
        ch_versions = ch_versions.mix(CIRCRNA_CROSSRUN_MERGE.out.versions.first())
    }

    //
    // circRNA quantification
    //

    if (runQuantify) {
        def ch_sample_discovery = CIRCRNA_ANALYSIS.out.discovery_bed
            .join(CIRCRNA_ANALYSIS.out.discovery_conf, by: 0)
            .map { meta, bed, tsv -> [meta, bed, tsv] }

        def ch_unit_catalog
        def ch_run_to_unit

        if (runCrossrunMerge) {
            def ch_group_discovery = CIRCRNA_CROSSRUN_MERGE.out.bed
                .filter { m, _bed -> m.tier == 'discovery' }
                .join   ( CIRCRNA_CROSSRUN_MERGE.out.confidence.filter { m, _c -> m.tier == 'discovery' }, by: 0 )

            ch_unit_catalog = ch_group_discovery
                .map { m, bed, conf -> [m.id, bed, conf] }

            def ch_sid_to_unit = ch_group_discovery
                .flatMap { m, _bed, _conf -> m.sample_ids.collect { sid -> [sid, m.id] } }

            ch_run_to_unit = ch_fastq
                .map { meta, _fq -> [meta.id, meta] }
                .join( ch_sid_to_unit, by: 0 )
                .map { _sid, meta, unit_id -> [meta, unit_id] }
        } else {
            ch_unit_catalog = ch_sample_discovery
                .map { meta, bed, conf -> [meta.id, bed, conf] }

            ch_run_to_unit = ch_fastq
                .map { meta, _fq -> [meta, meta.id] }
        }

        CIRCRNA_QUANTIFY (
            ch_fastq,
            ch_sample_discovery,
            ch_unit_catalog,
            ch_run_to_unit,
            file(params.fasta, checkIfExists: true),
            CIRCRNA_ANALYSIS.out.isocirc_native,
            CIRCRNA_ANALYSIS.out.circfl_native,
            CIRCRNA_ANALYSIS.out.cirilong_reads,
            CIRCRNA_ANALYSIS.out.circnick_native
        )
        ch_versions = ch_versions.mix(CIRCRNA_QUANTIFY.out.versions)

        def ch_clean_for_quant
        if (runCrossrunMerge) {
            ch_clean_for_quant = ch_run_to_unit
                .map     { meta, unit_id -> [unit_id, meta] }
                .combine ( CIRCRNA_CROSSRUN_MERGE.out.clean.map { m, tsv -> [m.id, m.tier, tsv] }, by: 0 )
                .map     { _unit_id, meta, tier, tsv -> [meta + [category: tier], tsv] }
        } else {
            ch_clean_for_quant = CIRCRNA_ANALYSIS.out.clean_tsv
                .filter { meta, _tsv -> meta.category in ['discovery', 'balanced_precision', 'balanced_recall', 'high_confidence'] }
        }

        def ch_clean_for_append = ch_clean_for_quant
            .map { meta, tsv -> [meta.id, meta, tsv] }
            .combine( CIRCRNA_QUANTIFY.out.cluster_map.map { meta, cmap -> [meta.id, cmap] }, by: 0 )
            .combine( CIRCRNA_QUANTIFY.out.final_counts.map { meta, counts -> [meta.id, counts] }, by: 0 )
            .map { _id, meta, tsv, cmap, counts -> [meta, tsv, cmap, counts] }

        QUANT_APPEND_COUNTS ( ch_clean_for_append )
        ch_versions = ch_versions.mix(QUANT_APPEND_COUNTS.out.versions.first())

        def ch_bed_for_quant
        if (runCrossrunMerge) {
            ch_bed_for_quant = ch_run_to_unit
                .map     { meta, unit_id -> [unit_id, meta] }
                .combine ( CIRCRNA_CROSSRUN_MERGE.out.bed.map { m, bed -> [m.id, m.tier, bed] }, by: 0 )
                .map     { _unit_id, meta, tier, bed -> [meta + [category: tier], bed] }
        } else {
            ch_bed_for_quant = CIRCRNA_ANALYSIS.out.discovery_bed.map        { m, b -> [m + [category: 'discovery'],          b] }
                .mix( CIRCRNA_ANALYSIS.out.balanced_bed.map        { m, b -> [m + [category: 'balanced_precision'],  b] } )
                .mix( CIRCRNA_ANALYSIS.out.balanced_recall_bed.map { m, b -> [m + [category: 'balanced_recall'],     b] } )
                .mix( CIRCRNA_ANALYSIS.out.high_conf_bed.map       { m, b -> [m + [category: 'high_confidence'],    b] } )
        }

        def ch_for_confident_filter = QUANT_APPEND_COUNTS.out.clean
            .filter { meta, _tsv -> meta.category in ['discovery', 'balanced_recall', 'high_confidence'] }
            .map    { meta, tsv -> [[meta.id, meta.category], meta, tsv] }
            .combine(
                ch_bed_for_quant
                    .filter { meta, _bed -> meta.category in ['discovery', 'balanced_recall', 'high_confidence'] }
                    .map    { meta, bed -> [[meta.id, meta.category], bed] },
                by: 0
            )
            .map { _key, meta, tsv, bed -> [meta, tsv, bed] }

        FILTER_CONFIDENT_DISCOVERY ( ch_for_confident_filter )
        ch_versions = ch_versions.mix(FILTER_CONFIDENT_DISCOVERY.out.versions.first())

        // balanced_precision never goes through FILTER_CONFIDENT_DISCOVERY, so
        // its final clean_with_counts.tsv is still QUANT_APPEND_COUNTS.out.clean
        def ch_final_clean_with_counts = QUANT_APPEND_COUNTS.out.clean
            .filter { meta, _tsv -> meta.category == 'balanced_precision' }
            .mix( FILTER_CONFIDENT_DISCOVERY.out.clean )

        def ch_deseq2_input = ch_final_clean_with_counts
            .map { meta, tsv -> [meta.category, meta, tsv] }
            .groupTuple(by: 0)

        BUILD_DESEQ2_MATRIX ( ch_deseq2_input )
        ch_versions = ch_versions.mix(BUILD_DESEQ2_MATRIX.out.versions)
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'nanocirc_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    if (!skipMultiqc) {
        def ch_multiqc_config        = channel.fromPath(
            "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        def ch_multiqc_custom_config = params.multiqc_config ?
            channel.fromPath(params.multiqc_config, checkIfExists: true) :
            channel.empty()
        def ch_multiqc_logo          = params.multiqc_logo ?
            channel.fromPath(params.multiqc_logo, checkIfExists: true) :
            channel.empty()

        def summary_params      = paramsSummaryMap(
            workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        def ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description                = channel.value(
            methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true
            )
        )

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
    }

    emit:
    multiqc_report = skipMultiqc ? [] : MULTIQC.out.report.toList()
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
