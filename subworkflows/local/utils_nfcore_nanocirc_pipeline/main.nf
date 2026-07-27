//
// Subworkflow with functionality specific to the nanocirc pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    def before_text = """
\033[0;32m                                         .. . .77?               .  . \033[0m
\033[0;32m .. .::.    :!:   ... :::    .:^~~^.  ^^:!?Y^.^?Y~  .7Y:!?Y!  ^^:!?Y. \033[0m
\033[0;32m:7^.::!?^   ^77.  :7:::^77:.^.   :77:5!    ..~JJ~.   5P..^P7:5~   ..  \033[0m
\033[0;32m^7     ~?  :^.77. !~     7!!:     :75B      .G!      BJ!P!  GP        \033[0m
\033[0;32m~?.    ^! .^:.:7! !!     !~!!:.  .:.7@?.  ..?@?.  ..:&5.~YJ:J&7.  .:  \033[0m
\033[0;32m.~~.   ~.::     7!:!^   .~ .^!!^:.   ~5Y7~:..JGPY?7^.YG.  ~555PJ!^:.  \033[0m
\033[0;32m      ^. .      ::     .^                                   ~5PY~.  ..\033[0m
\033[0;32m    .:.               .:                                      :?G#PJ~.\033[0m
\033[0;32m   ..                ..                                          ^:    \033[0m
-\033[2m----------------------------------------------------\033[0m-
\033[0;35m  nanocirc ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    def after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md
"""
    def command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq ->
                return [ meta + [ single_end:true ], fastq ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    def summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete = {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError = {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    return input
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citations = []
    citations << (params.skip_fastqc   ? "" : "FastQC (Andrews 2010)")
    citations << (params.skip_nanoplot ? "" : "NanoPlot (De Coster et al. 2018)")
    citations << (params.run_isocirc   ? "isoCirc (Xin et al. 2021)"        : "")
    citations << (params.run_circfl    ? "CircFL-seq (Liu et al. 2021)"     : "")
    citations << (params.run_cirilong  ? "CIRI-long (Zhang et al. 2021)"    : "")
    citations << (params.run_circnick  ? "circnick-lrs (Rahimi et al. 2021)" : "")
    citations << "bedtools (Quinlan & Hall 2010)"
    citations << (params.run_crossrun_merge ? "pybedtools (Dale et al. 2011)" : "")
    citations << (params.circnick_liftover_chain ? "UCSC liftOver (Kent et al. 2002)" : "")
    citations << (params.skip_annotation ? "" : "samtools (Danecek et al. 2021)")
    citations << (params.skip_annotation ? "" : "GffRead/GffCompare (Pertea & Pertea 2020)")
    citations << (params.skip_annotation ? "" : "AGAT (Dainat)")
    citations << (params.run_quantify ? "minimap2 (Li 2018)" : "")
    citations << (params.run_quantify ? "BLAT (Kent 2002)" : "")
    citations << (params.run_quantify ? "BWA (Li & Durbin 2009)" : "")
    citations << (params.skip_multiqc ? "" : "MultiQC (Ewels et al. 2016)")

    def citation_text = "Tools used in the workflow included: " +
        citations.findAll { it }.join(', ') + "."

    return citation_text
}

def toolBibliographyText() {
    def refs = []
    refs << (params.skip_fastqc   ? "" : "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>")
    refs << (params.skip_nanoplot ? "" : "<li>De Coster W, D'Hert S, Schultz DT, Cruts M, Van Broeckhoven C. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. 2018;34(15):2666-2669. doi: 10.1093/bioinformatics/bty149</li>")
    refs << (params.run_isocirc   ? "<li>Xin R, Gao Y, Gao Y, et al. isoCirc catalogs full-length circular RNA isoforms in human transcriptomes. Nat Commun. 2021;12(1):266. doi: 10.1038/s41467-020-20459-8</li>" : "")
    refs << (params.run_circfl    ? "<li>Liu Z, Tao C, Li S, et al. circFL-seq reveals full-length circular RNAs with rolling circular reverse transcription and nanopore sequencing. eLife. 2021;10:e69457. doi: 10.7554/eLife.69457</li>" : "")
    refs << (params.run_cirilong  ? "<li>Zhang J, Hou L, Zuo Z, et al. Comprehensive profiling of circular RNAs with nanopore sequencing and CIRI-long. Nat Biotechnol. 2021;39(7):836-845. doi: 10.1038/s41587-021-00842-6</li>" : "")
    refs << (params.run_circnick  ? "<li>Rahimi K, Venø MT, Dupont DM, Kjems J. Nanopore sequencing of brain-derived full-length circRNAs reveals circRNA-specific exon usage, intron retention and microexons. Nat Commun. 2021;12(1):4825. doi: 10.1038/s41467-021-24975-z</li>" : "")
    refs << "<li>Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010;26(6):841-842. doi: 10.1093/bioinformatics/btq033</li>"
    refs << (params.run_crossrun_merge ? "<li>Dale RK, Pedersen BS, Quinlan AR. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. Bioinformatics. 2011;27(24):3423-3424. doi: 10.1093/bioinformatics/btr539</li>" : "")
    refs << (params.circnick_liftover_chain ? "<li>Kent WJ, Sugnet CW, Furey TS, et al. The human genome browser at UCSC. Genome Res. 2002;12(6):996-1006. doi: 10.1101/gr.229102</li>" : "")
    refs << (params.skip_annotation ? "" : "<li>Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. Gigascience. 2021;10(2):giab008. doi: 10.1093/gigascience/giab008</li>")
    refs << (params.skip_annotation ? "" : "<li>Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. F1000Res. 2020;9:ISCB Comm J-304. doi: 10.12688/f1000research.23297.2</li>")
    refs << (params.skip_annotation ? "" : "<li>Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. Zenodo. doi: 10.5281/zenodo.3552717</li>")
    refs << (params.run_quantify ? "<li>Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191</li>" : "")
    refs << (params.run_quantify ? "<li>Kent WJ. BLAT-the BLAST-like alignment tool. Genome Res. 2002;12(4):656-664. doi: 10.1101/gr.229202</li>" : "")
    refs << (params.run_quantify ? "<li>Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics. 2009;25(14):1754-1760. doi: 10.1093/bioinformatics/btp324</li>" : "")
    refs << (params.skip_multiqc ? "" : "<li>Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048. doi: 10.1093/bioinformatics/btw354</li>")

    def reference_text = refs.findAll { it }.join(' ')

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
