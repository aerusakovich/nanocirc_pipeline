# nanocirc: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the top-level results directory specified with `--outdir`.

## Pipeline overview

The pipeline processes long-read nanopore FASTQ files through the following steps:

1. **Quality control**: FastQC and NanoPlot assess read quality
2. **circRNA detection**: up to four tools run in parallel (isoCirc, CircFL-seq, CIRI-long, circnick-lrs)
3. **BED12 conversion**: each tool's output is converted to a unified BED12 format
4. **Merging & confidence scoring**: when two or more tools are active, results are merged using the hybrid smart-merge algorithm and scored on two independent confidence axes
5. **Annotation**: each merged output is annotated with GFF comparison, FASTA extraction, circRNA type classification, and expression counts
6. **MultiQC**: aggregated QC report

---

## Quality control

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `qc/<sample>/fastqc/`
  - `*.html`: FastQC report with quality metrics per sample
  - `*.zip`: Archive containing the report and raw data

</details>

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports per-read quality scores, GC content, sequence length distribution, and adapter content. Useful for flagging low-quality or contaminated samples before analysis.

### NanoPlot

<details markdown="1">
<summary>Output files</summary>

- `qc/<sample>/nanoplot/`
  - `NanoPlot-report.html`: Interactive HTML report
  - `NanoStats.txt`: Summary statistics (N50, mean length, mean quality, etc.)
  - Various PNG plots: read length histogram, quality vs length scatter, etc.

</details>

[NanoPlot](https://github.com/wdecoster/NanoPlot) is designed specifically for nanopore data and provides read length distributions, quality score distributions, and yield-over-time plots.

---

## circRNA detection

Each active tool writes its raw output to a dedicated subdirectory and its BED12 file to a shared `bed12/` directory.

### isoCirc

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/isocirc/`
  - `isocirc_output/isocirc.bed`: Main isoCirc output (BED format with isoform detail)
  - `isocirc_output/`: Full isoCirc output directory

</details>

[isoCirc](https://github.com/Xinglab/isoCirc) detects circRNA isoforms from nanopore long reads using full-length read alignment.

### CircFL-seq

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/circfl_seq/`
  - `circFL_final.bed`: Final CircFL-seq output
  - Full pipeline output directory including intermediate RG, DNSC, cRG, mRG steps

</details>

[CircFL-seq](https://github.com/yangence/circfull) reconstructs full-length circRNA isoforms from nanopore sequencing using a rolling-circle amplification model.

### CIRI-long

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/ciri_long/`
  - `<sample>.info`: Main CIRI-long output with circRNA calls and isoform structure
  - `<sample>.isoforms`: Isoform-level output

</details>

[CIRI-long](https://github.com/bioinfo-biols/CIRI-long) uses a seed-and-extend strategy for circRNA detection and isoform characterisation from long reads.

### circnick-lrs

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/circnick_lrs/`
  - `<sample>/<sample>.circRNA_candidates.annotated.txt`: Annotated circRNA candidates
  - `<sample>/<sample>.circ_circRNA_exon_usage_length_of_exons.txt`: Exon usage per circRNA
  - `<sample>/<sample>.introns...intronCov.bed`: Intron coverage file
- `circrna/<sample>/circnick_lrs/lifted/` _(only if `--circnick_liftover_chain` was provided)_
  - `*_lifted_annotated.txt`: Coordinates lifted to target genome build
  - `*_lifted_exon_usage.txt`
  - `*_lifted_intron_cov.bed`
  - `*_liftover_failed.tsv`: circRNAs excluded due to failed liftover

</details>

[circnick-lrs](https://github.com/dzhang32/circnick) uses built-in mm10 or hg19 references. Provide `--circnick_liftover_chain` if your analysis uses a different genome build.

### BED12 files

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/bed12/`
  - `<sample>_cirilong.bed12`: CIRI-long output in BED12
  - `<sample>_circnick.bed12`: circnick-lrs output in BED12
  - `isocirc.bed`: isoCirc output (already in BED format)
  - `circFL_final.bed`: CircFL-seq output (already in BED format)

</details>

All tool outputs are converted to a 12-column BED format (BED12) for downstream merging. BED12 columns: chrom, start, end, name, score (read count), strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts.

---

## Merged outputs

Merging is performed when **two or more** detection tools are active. All tools are first grouped by relaxed BSJ coordinates (within `--circrna_bsj_tolerance` bp). Within each group, a merge algorithm selects the representative BSJ and exon structure. The result is confidence-scored on two independent axes and optionally filtered before publication.

### Pairwise comparisons

All pairwise combinations of active tools are compared using `bedtools intersect -split -wo` to identify shared circRNA isoforms (used internally for isoform confidence scoring). These intermediate comparison files are not published to the results directory.

### Merge algorithms

All tools within a relaxed-BSJ group (coordinates within `--circrna_bsj_tolerance` bp) are treated as candidates for the same circRNA. A merge algorithm selects the representative BSJ and exon structure from those candidates.

**BSJ selection** (all modes except `priority`): majority vote across all tools; ties broken by tool priority IsoCirc > CircFL-seq > CircNick-LRS > CIRI-long.

**Structure selection** differs between modes and is described in the tables below.

#### Default merge mode

The pipeline uses `consensus_hybrid` as its merge algorithm.

| Property | `consensus_hybrid` |
| -------- | ------------------ |
| BSJ | Majority vote across all tools |
| Structure vote participants | Tools sharing the **exact winning BSJ** only |
| Coordinate comparison | Absolute genomic coords (boundaries within `--circrna_bsj_tolerance` bp) |
| Rebasing of minority-BSJ tools | **No**: tools with a different BSJ do not contribute to the structure vote |
| Minority-BSJ tool handling | Emitted as separate isoform entries at their own BSJ coordinates |
| Structure tie-break priority | IsoCirc > CircFL-seq > CIRI-long > CircNick-LRS |

This design selects the most-supported exon structure among tools that agree on the BSJ, without shifting coordinates from tools that landed at a slightly different junction. Minority-BSJ isoforms are preserved in the output rather than silently discarded.

#### Default output files

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/merged/`
  - `<sample>_discovery.bed12`: hybrid, unfiltered (maximum recall)
  - `<sample>_discovery_confidence.tsv`
  - `<sample>_balanced_precision.bed12`: hybrid + isocirc_only filter (F1 / precision-leaning - lowest FP for decent recall)
  - `<sample>_balanced_precision_confidence.tsv`
  - `<sample>_balanced_recall.bed12`: consensus algorithm + trusted_only filter (best F1, recall-leaning)
  - `<sample>_balanced_recall_confidence.tsv`
  - `<sample>_high_confidence.bed12`: hybrid + high_only_isocirc filter (best precision, trusted subset)
  - `<sample>_high_confidence_confidence.tsv`

</details>

Four confidence-filtered outputs are published by default. After merging, each entry is scored on two independent confidence axes (`bsj_consensus` and `isoform_consensus`, each Low / Medium / High) and one of the following filters is applied:

| Output | Merge algorithm | Filter | Rule | Axes retained |
| ------ | ---------------- | ------ | ---- | ------------- |
| **`discovery`** | `consensus_hybrid` | none | Keep all entries | any |
| **`balanced_precision`** | `consensus_hybrid` | `isocirc_only` | Drop entries where either axis is Low, unless the source is IsoCirc | ≥ Medium on both, or Low from IsoCirc |
| **`balanced_recall`** | `consensus` | `trusted_only` | Drop entries where either axis is Low, unless the source is CIRI-long or IsoCirc | ≥ Medium on both, or Low from a trusted tool |
| **`high_confidence`** | `consensus_hybrid` | `high_only_isocirc` | Keep only entries where both axes are High, unless Low and the source is IsoCirc | High on both, or Low from IsoCirc |

If `--run_quantify true` is also set, `discovery`/`balanced_recall`/`high_confidence` get a further post-quantification confidence filter on top of this table -- see [Quantification](#quantification) below. This filter does not run on `balanced_precision`: its "Low from IsoCirc" entries pass as-is, favoring recall, whereas `high_confidence` additionally requires read support there since that tier prioritizes precision over recall.

`balanced_recall` uses the `consensus` merge algorithm (string-equality structure comparison, see below) rather than `consensus_hybrid`. It is the same trusted_only-filtered output previously published only under `--run_benchmark_modes` as `smart_consensus_filtered`.

#### Additional merge modes (`--run_benchmark_modes`)

Two further algorithms remain available for research and benchmarking only. Both were outperformed by `consensus_hybrid`/`consensus` in benchmark evaluation and are not published by default.

| Mode | BSJ selection | Structure vote participants | Minority-BSJ rebasing | Minority-BSJ handling |
| ---- | ------------- | -------------------------- | --------------------- | --------------------- |
| `consensus_xstruct` | Majority vote | All tools in group | **Yes**, shifted to winning BSJ | Folded into structure vote |
| `priority` | Highest-priority tool | Highest-priority tool | No | Separate isoforms at own coords |

Key differences from `consensus_hybrid`:

- **`consensus`** (used by `balanced_recall` above) uses string equality for structure comparison rather than coordinate similarity: two tools reporting the same exon structure with a 1 bp boundary difference are treated as distinct isoforms.
- **`consensus_xstruct`** includes minority-BSJ tools in the structure vote by rebasing their exon coordinates to the winning BSJ. This can incorporate more structural information but may introduce coordinate imprecision when BSJ offset is non-trivial.
- **`priority`** skips voting entirely: BSJ and structure come unconditionally from the single highest-priority tool present.

With `--run_benchmark_modes`, `consensus` additionally publishes its `no_low` filter variant, and `consensus_xstruct`/`priority` each publish all four filter variants (unfiltered, `no_low`, `trusted_only`, `high_only`).

### Confidence TSV format

All `*_confidence.tsv` files share a common format. Confidence is assessed on two **independent axes**:

- **`bsj_consensus`**: what fraction of active tools agreed on this back-splice junction?
- **`isoform_consensus`**: how well do tools agree on the exon boundaries of this isoform?

Each axis is scored independently (1=Low, 2–3=Medium, 4=High based on binned percentage of tools) allowing a circRNA to have a well-supported BSJ but uncertain isoform structure, or vice versa.

| Column               | Description                                                                 |
| -------------------- | --------------------------------------------------------------------------- |
| `#chrom`             | Chromosome                                                                  |
| `start`              | BSJ start (0-based)                                                         |
| `end`                | BSJ end                                                                     |
| `strand`             | Strand (`+` or `-`)                                                         |
| `bsj_id`             | Unique identifier: `chrom:start-end:strand` (isoforms suffixed `\|iso*`)   |
| `bsj_confidence`     | Number of tools detecting this BSJ (1–4)                                    |
| `<tool>`             | Per-tool presence flag: `1` if detected, `0` if not (one column per tool)  |
| `<tool>_block_sizes` | BED12 block sizes from this tool's call                                     |
| `<tool>_block_starts`| BED12 block starts from this tool's call                                    |
| `isoform_confidence` | Number of tools with confirmed isoform overlap                              |
| `bsj_score`          | Percentage of active tools detecting this BSJ, binned 1–4                  |
| `isoform_score`      | Percentage of active tools with isoform support, binned 1–4 (min 1)        |
| `overlap_score`      | Average pairwise spliced-length overlap fraction, binned 1–4               |
| `bsj_consensus`      | BSJ confidence label: `Low` (score 1), `Medium` (2–3), `High` (4)         |
| `isoform_consensus`  | Isoform confidence label: `Low` (score 1), `Medium` (2–3), `High` (4)     |

**Scoring bins (percentage of active tools):**

| % of active tools | Score | Consensus |
| ----------------- | ----- | --------- |
| ≤ 25%             | 1     | Low       |
| ≤ 50%             | 2     | Medium    |
| ≤ 75%             | 3     | Medium    |
| > 75%             | 4     | High      |

> [!NOTE]
> Consensus labels always reflect agreement among the tools that were actually run.
> A `High` from 2 tools means both tools agreed. It is not mathematically equivalent
> to `High` from 4 tools. The pipeline emits a warning when fewer than 4 tools are active.

---

## Annotation

Each merged output (discovery, balanced_precision, balanced_recall, high_confidence) is annotated per sample. Annotation runs when `--skip_annotation` is not set (default: enabled).

### Output structure

<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/merged/`
  - `annotated/`: GFF comparison output and annotated TSV (with class codes and reference gene IDs)
    - `<sample>_<tier>.annotated.tsv`: Full annotated TSV including isoform rows
  - `clean/`: Wet-lab-friendly TSV (narrowed column set, type + expression added; same rows, including isoforms, as the annotated TSV)
    - `<sample>_<tier>_clean.tsv`
  - `fasta/`: FASTA sequences of detected circRNAs
    - `<sample>_<tier>.fa`
  - `gff/`: GFF3 files
    - `<sample>_<tier>.gff3`: GFF3 derived from BED12
    - `<sample>_<tier>.annotated.gtf`: GFFcompare annotated GTF

</details>

`<tier>` is one of `discovery`, `balanced_precision`, `balanced_recall`, or `high_confidence`.

### Clean TSV format

The clean TSV (`*_clean.tsv`) is the primary output for downstream analysis. It contains the same rows as the annotated TSV, including isoform rows (entries with `bsj_id` containing `|iso`, which represent minority-BSJ calls from tools that disagreed on the back-splice junction), with a narrowed column set plus `type`/`supporting_reads` added. `bsj_id` is never stripped of its isoform suffix, so it matches the annotated TSV and BED12 exactly.

| Column               | Description                                                                                      |
| -------------------- | ------------------------------------------------------------------------------------------------ |
| `#chrom`             | Chromosome                                                                                       |
| `start`              | BSJ start (0-based)                                                                              |
| `end`                | BSJ end                                                                                          |
| `strand`             | Strand (`+` or `-`)                                                                              |
| `sel_block_count`    | Number of exon blocks in the selected isoform                                                    |
| `sel_block_sizes`    | Comma-separated exon block sizes (bp)                                                            |
| `sel_block_starts`   | Comma-separated exon block starts relative to `start`                                           |
| `bsj_id`             | Unique identifier: `chrom:start-end:strand`                                                      |
| `bsj_confidence`     | Number of tools detecting this BSJ (1–4)                                                         |
| `isoform_confidence` | Number of tools with confirmed isoform overlap                                                   |
| `class_code`         | GFFcompare class code describing the relationship to the reference annotation                    |
| `ref_gene_id`        | Reference gene ID from the GTF (`.` if intergenic)                                              |
| `type`               | circRNA biotype: `eciRNA`, `EIciRNA`, `ciRNA`, `antisense`, or `intergenic` (see below)         |
| `supporting_reads`   | Read count from the highest-priority active tool (priority: isoCirc > CIRI-long > CircFL-seq > circnick-lrs) |
| `supporting_tools`   | Comma-joined names of every tool that called this locus (e.g. `isocirc,cirilong`). For a crossrun clean TSV, this is the union across every contributing sample's matched row. |

#### circRNA type classification

Types are assigned by intersecting BSJ coordinates against gene and exon BED files derived from the GTF:

| Type          | Definition                                                                          |
| ------------- | ----------------------------------------------------------------------------------- |
| `eciRNA`      | Same-strand gene overlap; the circRNA is fully covered by exonic regions (purely exonic) |
| `EIciRNA`     | Same-strand gene overlap; overlaps exons but retains intronic content (exon–intron circRNA) |
| `ciRNA`       | Same-strand gene overlap; no exon overlap (purely intronic)                         |
| `antisense`   | Overlaps a gene on the opposite strand only                                         |
| `intergenic`  | No overlap with any annotated gene on either strand                                 |

---

## Cross-run merge

When `--run_crossrun_merge true` is set and the samplesheet contains a `group` column, all samples sharing the same group are merged together after per-sample analysis. Each run is treated as an independent caller. The same `consensus_hybrid` algorithm and confidence scoring used within a sample (across tools) is applied across runs (across samples).

In plain terms, the four tiers mean:

| Tier              | What it takes to be retained |
| ----------------- | ----------------------------- |
| **`discovery`**   | Detected by **at least 1 tool** in **at least 1 run**: maximum sensitivity, use for exploration, understanding it has the highest FP rate of all modes |
| **`balanced_precision`** | **Multiple tools agreed** within a run AND **multiple runs** support it: recommended for precision leaning analyses |
| **`balanced_recall`** | **Multiple tools agreed** within a run AND **multiple runs** support it: recommended for recall leaning analyses |
| **`high_confidence`** | **All tools agreed** within runs AND **most/all runs** support it: maximum precision, lowest false positive rate |

### How cross-run merge decides a locus

Each candidate locus is scored with two separate votes, run independently:

1. **BSJ vote**: counts how many different runs support the exact back-splice junction position, regardless of exon structure. This vote decides which tier a locus qualifies for (the count thresholds above).
2. **Structure vote**: restricted to records at the winning BSJ position only, groups them by exon structure and ranks groups by total tool agreement. This vote decides the final isoform structure reported for the locus.

Splitting these two votes avoids a locus's correct BSJ call being discarded just because a different run happened to report a slightly different exon structure at the same junction.

A structure seen in only 1 run is dropped from the structure vote unless that run's own tool agreement meets `--crossrun_min_tool_agreement` (default `2`); structures backed by 2 or more runs are never dropped by this check. This removes weak, single-run structure calls without changing which tier the underlying BSJ qualifies for. Minority BSJ positions (backed by fewer runs than the winning position) each keep their own separate isoform entry rather than being merged away.

### Output files

<details markdown="1">
<summary>Output files</summary>

- `circrna/crossrun/<group>/<tier>/`
  - `<group>_<tier>_crossrun.bed12`: Merged BED12 filtered to the tier threshold
  - `<group>_<tier>_crossrun_confidence.tsv`: Full intermediate TSV: all merged circRNAs with per-sample BSJ and isoform consensus columns
  - `<group>_<tier>_crossrun_clean.tsv`: Wet-lab-friendly TSV, same core columns as the [Clean TSV format](#clean-tsv-format) (`#chrom` through `type`, `supporting_tools`) but without `class_code`/`ref_gene_id`/`supporting_reads`; adds `n_samples` and a `<run>_bsj_consensus`/`<run>_isoform_consensus` pair per contributing run

</details>

`<tier>` is one of `discovery`, `balanced_precision`, `balanced_recall`, or `high_confidence`. `<group>` is the group name from the samplesheet.

### Count thresholds

Minimum number of runs that must detect a circRNA for it to be retained, where `n` is the total number of runs in the group:

| Tier              | Minimum runs required        |
| ----------------- | ---------------------------- |
| `discovery`       | ≥ 1 (all circRNAs retained)  |
| `balanced_precision` | ≥ max(2, ceil(0.25 × n))  |
| `balanced_recall` | ≥ max(2, ceil(0.25 × n))     |
| `high_confidence` | ≥ ceil(0.75 × n)             |

### Cross-run confidence TSV format

The `*_crossrun_confidence.tsv` has the same core columns as the per-sample confidence TSV, plus:

| Column                        | Description                                              |
| ----------------------------- | -------------------------------------------------------- |
| `n_samples`                   | Number of runs detecting this circRNA                    |
| `type`                        | circRNA biotype (see [type classification](#circrna-type-classification)) |
| `<run>_bsj_consensus`         | BSJ consensus label for this run (`Low`/`Medium`/`High`, or empty if not detected) |
| `<run>_isoform_consensus`     | Isoform consensus label for this run                     |

---

## Quantification

When `--run_quantify true` is set, each sample's reads are remapped against a set of synthetic tandem-duplicated circle references built from a `discovery` (unfiltered, hybrid-consensus) locus catalog: a chunked remap-and-classify pass (tier1), an overlap-cluster rescue pass, then targeted low-coverage rescue (tier2) and gene-family/repeat-cluster rescue (tier3) for loci tier1 handled poorly.

**Which catalog a sample is quantified against depends on `--run_crossrun_merge`:**

- **Off**: each sample is quantified against its own per-sample `discovery` catalog.
- **On**: every sample in the same samplesheet `group` is quantified against **one shared catalog** built from that group's crossrun `discovery` output, so per-locus counts line up across samples in the group (needed for any downstream count-matrix comparison, e.g. differential expression). The circle references, chunk references, and gene-family similarity database are built once per group, not once per sample; only the actual read alignment/classification is per sample.

Quantification counts are also appended onto each sample's per-tier clean TSV (`circrna_clean.py`/`crossrun_annotate.py` output) as three new columns, producing a `_clean_with_counts.tsv` table:

| Column                           | Description                                                              |
| -------------------------------- | ------------------------------------------------------------------------ |
| `nanocirc_quant_reads`           | Final read count for this locus (0 if quantification found no count)     |
| `nanocirc_quant_tier`            | Which stage produced the count: `tier1`, `tier2`, `tier3`, or `tier1_gene_family_unresolved` |
| `nanocirc_quant_low_confidence`  | `true` for gene-family/repeat-cluster loci tier3 could not resolve further |

- `discovery` / `balanced_recall`: drops loci only CircNick-LRS called (`supporting_tools` is just `circnick`), with `nanocirc_quant_reads` at or below `--circrna_confident_min_reads`.
- `balanced_precision` normally requires at least 2 tools to agree, but it keeps Isocirc only calls, based on the strength of IsoCirc's high precision.
- `high_confidence`: normally requires most active tools to agree but also lets IsoCirc-only calls that are meeting read support threshold. This filter drops the IsoCirc-only entries with `nanocirc_quant_reads` at or below `--circrna_confident_min_reads`.



<details markdown="1">
<summary>Output files</summary>

- `circrna/<sample>/quantify/`
  - `<sample>_quant_final_counts.tsv`: raw per-locus read counts from the tier1/tier2/tier3 rescue passes, before being joined onto the clean TSV
- `circrna/<sample>/merged/clean/`
  - `<sample>_<tier>_clean_with_counts.tsv`: clean TSV with quantification columns appended (see above), and (for `discovery`/`balanced_recall`/`high_confidence`) the confidence filter applied
- `circrna/<sample>/merged/`
  - `<sample>_<tier>.bed12`: matching BED12, same filtering applied where relevant (supersedes the pre-quantification version at this same path)

</details>

All other intermediate quantification files (reference chunks, genome-wide alignments, the similarity database, per-tier raw counts before the final merge) are internal to each task's work directory and are not published.

Quantification-specific parameters:

| Parameter                       | Description                                                              | Default |
| -------------------------------- | ------------------------------------------------------------------------ | ------- |
| `--run_quantify`                 | Enable quantification                                                     | `false` |
| `--quant_chunk_size`             | Loci per reference chunk for tier1                                       | `100`   |
| `--quant_chunk_seed`             | Fixed shuffle seed for chunk assignment                                  | `42`    |
| `--quant_locus_dedup_tolerance`  | bp tolerance (both start AND end) for near-duplicate locus dedup         | `10`    |
| `--quant_min_old_tool_count`     | Min independent old-tool-count to flag a locus for tier2 rescue          | `50`    |
| `--circrna_confident_min_reads`  | Max quantified read count for the discovery/balanced_recall confidence filter (see above) | `2` |

### DESeq2 count matrix

When `--run_quantify true` is set, the pipeline also builds one wide isoform x sample count matrix per confidence tier, pooling every sample in the run. A feature row is one isoform: two isoforms sharing a BSJ but with different exon structure get separate rows. Samples quantified against the same crossrun catalog (same `group`, `--run_crossrun_merge true`) already share row identity; samples from different catalogs are unioned, 0-filled where a sample's own catalog never called that isoform.

<details markdown="1">
<summary>Output files</summary>

- `circrna/deseq2/`
  - `deseq2_counts_<tier>.tsv`: isoform x sample raw read counts (`nanocirc_quant_reads`), ready for `DESeqDataSetFromMatrix`
  - `deseq2_coldata_<tier>.tsv`: `sample`, `group` columns, same row order as the counts matrix's columns
  - `deseq2_features_<tier>.tsv`: per-isoform coordinates, exon block structure, `bsj_id`, `type`

</details>

---

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: Standalone HTML report viewable in any browser
  - `multiqc_data/`: Parsed statistics from all tools
  - `multiqc_plots/`: Static plot images

</details>

[MultiQC](https://multiqc.info) aggregates QC results from FastQC and NanoPlot across all samples into a single report. Skip with `--skip_multiqc`.

---

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_*.html`: Nextflow execution report (resource usage per process)
  - `execution_timeline_*.html`: Timeline of all processes
  - `execution_trace_*.txt`: Raw trace file with per-task metrics
  - `pipeline_dag_*.html`: Directed acyclic graph of the pipeline
  - `nf_core_nanocirc_software_mqc_versions.yml`: Software versions for all tools

</details>

Nextflow automatically generates execution reports for every run. These are useful for troubleshooting, optimising resource requests, and recording the exact software versions used.
