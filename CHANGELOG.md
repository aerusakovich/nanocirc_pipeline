# nf-core/nanocirc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/nanocirc, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- New optional quantification stage (`--run_quantify`): remap-based circRNA read counting. Without `--run_crossrun_merge`, each sample is quantified against its own discovery catalog. With it, every sample in a `group` shares one catalog built from that group's crossrun discovery output. New container `quay.io/anrusakovich/nanocirc-quant`.
- New confidence filter (`FILTER_CONFIDENT_DISCOVERY`, `bin/filter_confident_discovery.py`), applied to the `discovery` and `balanced_recall` tiers after quantification to remove a share of FP. Drops two low-confidence patterns with weak read support: loci only CircNick-LRS called, and intergenic/antisense loci called by only one tool if they didn't have enough read support. New param `--circrna_confident_min_reads` (default `2`) sets the read-count cutoff. `balanced_precision`/`high_confidence` are untouched (their own stricter multi-tool rules already exclude this pattern). Part of the new quantify. 
- New `supporting_tools` column in clean TSVs (`bin/circrna_clean.py`, `bin/crossrun_annotate.py`): comma-joined names of every tool that called a locus (union across samples for cross-run merge). Needed by the new confidence filter above, and useful on its own for downstream filtering.
- Redesigned cross-run merge (`bin/smart_merge.py`'s `cross_run_hybrid_entries()`) to score each locus with two separate votes instead of one: a BSJ-position vote (decides which tier a locus qualifies for) and a structure vote restricted to the winning BSJ position (decides the final isoform structure). This stops a locus's correct BSJ call from being discarded just because a different run reported a slightly different exon structure at the same junction. Minority BSJ positions now keep their own separate isoform entry instead of being merged away. New param `--circrna_crossrun_min_corroboration` (default `2`) drops a single-run structure call unless that run's own tool agreement meets this bar; this cut cross-run false positives by about 34% in testing, for a small recall cost.
- Promoted the trusted_only-filtered `consensus` merge output out of `--run_benchmark_modes` into the default pipeline as a permanent 4th confidence tier, `balanced_recall` (previously published only as `smart_consensus_filtered` in benchmark mode).

### `Changed`

- Renamed the `balanced` confidence tier to `balanced_precision` for clarity now that a recall-leaning counterpart (`balanced_recall`) exists alongside it (per-sample outputs, cross-run merge tiers, and docs updated accordingly).
- Updated `bin/smart_merge.py`'s `BSJ_PRIORITY`/`STRUCT_PRIORITY` tie-break orders to reflect updated benchmark results.
- Flattened merged output layout: `merged/smart/` and `merged/pairs/` no longer exist as separate subdirectories. Smart-merge/filter/annotate/finalize outputs publish directly under `merged/` (legacy `strict`/`relaxed`/`exon_based` benchmark-mode subdirectories are unchanged); pairwise bedtools comparison files are no longer published.

### `Fixed`

- `QUANT_GENOME_INDEX` bumped from `process_low` (12GB) to `process_high` (72GB, `* task.attempt`) with explicit retry, confirmed by a real SIGKILL (OOM) failure building a minimap2 splice index for a full mouse genome.
- Boolean CLI flags (`--run_isocirc`, `--run_circfl`, `--run_cirilong`, `--run_circnick`, `--run_crossrun_merge`, `--run_quantify`, `--run_benchmark_modes`, `--skip_qc`, `--skip_fastqc`, `--skip_nanoplot`, `--skip_multiqc`, `--skip_annotation`) passed with an explicit `false` (e.g. `--run_circnick false`) had no effect: `NXF_SYNTAX_PARSER=v2` delivers CLI params as strings, and a non-empty string is truthy in Groovy, so every `if (params.run_x)` gate ran regardless of what was passed. A previous fix attempt (`params.x = params.x.toBoolean()` in `main.nf`) did not work either, confirmed with a minimal repro: Nextflow keeps returning the original CLI string on every later read of `params.x`, ignoring the reassignment. Found via a smoke test that toggled tools off one at a time: disabling `circnick` still ran `CIRCNICK_LRS`, and toggling off `run_crossrun_merge` still ran `CIRCRNA_CROSSRUN_MERGE`, which then crashed on a single-sample group (`add_isoform_confidence.py --pairs` needs at least 2). Fixed by reading each flag once into a local boolean (`asBool()` in `workflows/nanocirc.nf` and `subworkflows/local/circrna_analysis.nf`, inlined in `conf/modules.config`) instead of relying on `params.x`'s truthiness directly.
- `CIRCRNA_FINALIZE` crashed with an "input file name collision" whenever exactly 2 of the 4 detection tools were inactive at once: their expression-file placeholders (`ch_iso_expr`/`ch_fl_expr`/`ch_nick_expr`/`ch_ciri_expr` in `subworkflows/local/circrna_analysis.nf`) all used the same literal filename `NO_FILE`, and Nextflow can't stage two different missing-input placeholders under an identical name in the same task. Only surfaced once the boolean-flag fix above made turning tools off actually work. Gave each tool its own placeholder name (`NO_FILE_ISOCIRC`/`NO_FILE_CIRCFL`/`NO_FILE_CIRCNICK`/`NO_FILE_CIRILONG`, matching the convention already used elsewhere in the same file and in `quant_coverage_flag.nf`), and updated `modules/local/circrna_finalize.nf`'s presence checks to match.

### `Dependencies`

### `Deprecated`
