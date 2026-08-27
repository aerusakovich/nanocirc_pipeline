# nanocirc: Benchmark Results and Recommendations

This benchmark is adapted from our own earlier benchmark study, Rusakovich et al., "Benchmarking circRNA detection tools from long-read sequencing using a data-driven and flexible simulation framework" (Peer Community Journal, 2026; see [CITATIONS.md](../CITATIONS.md) for the full reference), modified and extended based on our newer findings. 

This page summarizes the benchmark behind nanocirc's three confidence tiers and gives practical guidance on which one to use. The new simulator and benchmark scripts is in this repo, [nano-circ](https://gitlab.com/bingo-igdr/nano-circ), as v3. Results below are from simulated human and mouse datasets, 3 runs each.

Several benchmarks are shown below:

- **[Same circRNAs, three wet-lab protocols](#same-circrnas-three-wet-lab-protocols)**: NanoSim's read-error model is characterized from real CIRI-long/isoCIRC and circFL-seq protocols ONT control reads; circRNA type composition and length/exon/splice distributions are the same across all 3 wet-lab methods, circRNA generation based on circRNAs from a human circbase+circatlas catalog. 
- **[Species benchmark: human vs mouse](#species-benchmark-human-vs-mouse)**: CIRI-long protocol, each species' own circRNA catalog and genome. Any difference here comes from species biology (genome, catalog composition), not the wet-lab protocol, which is fixed to CIRI-long for both.
- **[Database benchmark: annotated circRNA recovery](#database-benchmark-annotated-circrna-recovery)**: annotated vs de novo circRNA recovery. Same simulation and pipeline as above, but circRNAs are sampled from a real human circBase/circAtlas catalog instead of generated de novo, so tools are evaluated on already-known, already-catalogued circRNAs instead of novel unannotated ones.

## Recommendations

- **Choose circFL-seq wet-lab protocol if you can, followed by CIRI-long, for maximum tool compatibility.** Both work with all 4 tools; circFL-seq protocol gave the highest recall, highest F1, and the highest BSJ-confirmation rate of the 3 protocols we tested (59.2% vs CIRI-long's 51.7%), and its read structure also allowed to better resolve isoform structure. isoCirc-protocol reads are the hardest for the tools to work with, and CIRI-long specifically is almost incompatible with them. See [Structural precision/recall](#structural-precisionrecall) and [Isoform detection](#isoform-detection) under [Same circRNAs, three wet-lab protocols](#same-circrnas-three-wet-lab-protocols).
- **Use nanocirc itself rather than any single tool, and use its own quantification (`--run_quantify`) rather than a raw tool's built-in read counting.** All 4 tools are containerized and patched against real bugs found running them at scale, several silent rather than crashing (wrong fastq/genome argument that always caused crash, a BAM filename mismatch that silently dropped intron-coverage output, a BioPython-version crash that only crashes on certain read alignments). A tool being good at finding the right circRNA doesn't mean its own built-in quantifier also counts it well and native quantifiers can't quantify circRNAs tool didn't report. nanocirc's remap-based quantification is not dependent on which tool/tier found a locus. It applies the same counting method regardless of source, is stable across wet-lab approaches and has been benchmarked to agree well across tools/tiers and against ground truth (see [Quantification advice](#quantification-advice)). It also builds one shared, cross-sample circRNA catalog (`--run_crossrun_merge`) for a multi-sample count matrix that can be used with DESeq2 directly.
- **For isoform-level differential expression, use the `balanced` tier.** It is close to `discovery`'s isoform-family resolution while keeping meaningfully higher precision, so you get `discovery`'s isoform sensitivity with fewer false positives in the count matrix. `discovery` is an acceptable choice for locus-level differential expression, where getting the exact internal exon structure right isn't the point. See [Isoform detection](#isoform-detection).
- **Exception: on isoCirc-protocol data, use `discovery` for isoform-level work instead of `balanced`.** circNICK-lrs's own structure calls aren't trustworthy enough to count as consensus evidence, and CIRI-long is close to incompatible with isoCirc-protocol reads (see [Which tools work with my wet-lab protocol and species?](#which-tools-work-with-my-wet-lab-protocol-and-species)), so the only tools left to agree on a recovered isoform are circFL-seq and isoCirc, which rarely call the exact same structure independently. See [Isoform detection](#isoform-detection).
- **To verify a specific circRNA's exon boundaries before follow-up work (primers, probes, isoform-specific validation), check it against the `high_confidence` set.** It only contains calls with strong cross-tool agreement on both BSJ and structure, so an entry present there has cleared more quality thresholds. See [Which tier should I use?](#which-tier-should-i-use) and [Isoform detection](#isoform-detection).

## Which tools work with my wet-lab protocol and species?

nanocirc runs 4 tools (CIRI-long, circFL-seq, isoCirc, circNICK-lrs) by default. Not every tool works with every wet-lab protocol's read structure. Set `--wet_lab` to your protocol (`ciri_long`, `circfl`, `isocirc`, or `circnick`) to use default settings for this type of data. `--wet_lab circnick` switches to circNICK-lrs alone: this is the one protocol where the other tools would fail, since they expect rolling-circle reads with multiple copies of the circRNA per read, which is not the case for circNICK-lrs wet-lab data. `--wet_lab isocirc` turns CIRI-long off by default (see below). For every other protocol, `--wet_lab` does not turn any tool off by itself: so we suggest you to try with the default settings, but you can turn on/off tools with `--run_isocirc`, `--run_circfl`, `--run_circnick`, `--run_cirilong` set to `false`.

| Your wet-lab protocol | Tools suggestion |
| --- | --- |
| circNICK-lrs | circNICK-lrs only |
| CIRI-long | All 4 tools |
| circFL-seq | All 4 tools |
| isoCirc* | 3 tools, CIRI-long off by default** |
| PacBio | protocol-dependent, turn CIRI-long off manually if isoCirc based*** |
| Species | protocol-dependent, turn off circNICK-lrs for any species other than human/mouse**** |

\* circFL-seq also needs more memory than usual on isoCirc data: a first attempt at 250GB ran out of memory, a second attempt completed at 500GB, after about 5.5 days. Give circFL-seq extra memory if your samples come from isoCirc-protocol reads.

\*\* On real isoCirc-protocol CIRI-long did not finish in reasonable time, getting stuck on its collapse step for days. We tried it 3 times (5 days, then 10 days, then over 24 days), it still has not completed. It worked on simulated data and 10 pct of the real isoCirc dataset, but much slower than on the datasets from other wet-lab protocols, so this issue may be dataset size dependent. To try CIRI-long on isoCirc data use `--run_cirilong true`. By our benchmark, leaving it off on isoCirc data costs at most 0.01 F1 at any of the 4 tiers, while improving compute time and reducing the risk of pipeline freeze on this tool. 

\*\*\* On real isoCirc-based PacBio human data, CIRI-long did not finish in reasonable time, getting stuck on its collapse step for days. We tried it 3 times (5 days, then 10 days, then over 24 days), it still has not completed. It worked on simulated data and 10 pct of the real isoCirc-based PacBio dataset, so this issue may be dataset size dependent. PacBio has no dedicated `--wet_lab` preset, as it is a sequencing platform, so select the preset based on the wet-lab approach and turn `--run_cirilong` off manually if you hit the same issue.

\*\*\*\* circNICK-lrs is limited by species: it takes a hard-coded `--circnick_species` value of `human` or `mouse` only, and will not run for any other species, by design, whatever your wet-lab protocol is. Turn it off with `--run_circnick false` if you're not working with human or mouse data. If your protocol is circNICK-lrs's own and your species is neither human nor mouse, no tool in nanocirc will work for your data, since circNICK-lrs is the only tool compatible with that protocol.

## Which tier should I use?

Each tier balances recall against precision differently. There is no single best tier, only the right one for your question.

| Tier | What it optimizes for | Trade-off |
| ---- | ---------------------- | --------- |
| **`discovery`** | Recall. Keeps almost every circRNA any tool found. | Highest false-positive rate: computational artefacts and incorrectly predicted structures will pass alongside real calls. |
| **`balanced`** | highest F1, leaning recall. Removes most of `discovery`'s false-positive rate while keeping reasonable recall, including on isoform structure (see [Isoform detection](#isoform-detection)). | Good default for exploratory analyses and for isoform-level questions, where `discovery` alone would carry more false positives for similar isoform resolution. |
| **`high_confidence`** | Precision. Only calls with strong cross-tool agreement on both BSJ and structure, or isoCirc calls that met a defined read threshold, pass. | Trades some recall for precision compared to `balanced`, while keeping reasonable recall. Best for analyses leaning toward precision, not a maximally strict confirmatory set. |

> **Recommendation:** use **`balanced`** for most applications, since it has reached the highest F1. Reserve `discovery` for exploratory studies, and `high_confidence` for confirmatory ones that need a small, highly precise set.

## Quantification advice

Quantification is made by our pipeline's remap step, independent of which tool originally called a locus: the reproducibility and expression correlation plots below show read counts agree well across tiers, and against ground truth on remaining stable across all tested conditions.

For a multi-sample design, turn on `--run_crossrun_merge` and `--run_quantify` so every sample is evaluated within the `group` and is quantified against the shared group circRNA catalog. It allows to detect and quantify circRNAs that tools may have missed in one run despite the read evidence, but found in its replicates. Otherwise each sample will get its own per-run catalog.

## Same circRNAs, three wet-lab protocols

We used one circRNA catalog (same set of simulated circRNAs, human genome, same abundances) and only changed the wet-lab protocol NanoSim is trained on. Any difference in the results below comes from the wet-lab protocol characteristics.

NanoSim for each protocol was trained on:
CIRI-long protocol: Karpas ALCL (SRR18213998, Fuchs et al. 2022);
circFL-seq protocol: HEK293 (SRR14263341, Liu et al. 2021);
isoCirc protocol: HEK293 (Xin et al. 2021).

<details>
<summary><strong>Structural precision/recall</strong></summary>

### Structural precision/recall

<img src="images/benchmark/pr_scatter_equal_ciri_long.png" width="100%"/>
<img src="images/benchmark/pr_scatter_equal_circfl.png" width="100%"/>
<img src="images/benchmark/pr_scatter_equal_isocirc.png" width="100%"/>

Same circRNA catalog, same callers, only the protocol changed between these 3 plots.

Mean precision (P) / recall (R) / F1 across 3 runs, exon-based matching:

| Tier / tool | CIRI-long protocol | circFL-seq protocol | isoCirc protocol (CIRI-long on) |
| --- | --- | --- | --- |
| isoCirc alone | P 0.86 / R 0.14 / F1 0.24 | P 0.83 / R 0.18 / F1 0.30 | P 0.84 / R 0.08 / F1 0.14 |
| circFL-seq alone | P 0.82 / R 0.21 / F1 0.34 | P 0.80 / R 0.25 / F1 0.39 | P 0.54 / R 0.18 / F1 0.27 |
| CIRI-long alone | P 0.73 / R 0.33 / F1 0.46 | P 0.76 / R 0.41 / F1 0.53 | P 0.74 / R 0.03 / F1 0.05 |
| circNICK-lrs alone | P 0.19 / R 0.07 / F1 0.10 | P 0.19 / R 0.08 / F1 0.11 | P 0.17 / R 0.08 / F1 0.11 |
| `discovery` | P 0.67 / R 0.41 / F1 0.51 | P 0.67 / R 0.50 / F1 0.57 | P 0.57 / R 0.23 / F1 0.33 |
| `balanced` | P 0.76 / R 0.41 / F1 0.53 | P 0.77 / R 0.50 / F1 0.60 | P 0.71 / R 0.22 / F1 0.34 |
| `high_confidence` | P 0.86 / R 0.24 / F1 0.37 | P 0.86 / R 0.29 / F1 0.43 | P 0.82 / R 0.07 / F1 0.13 |

</details>

<details>
<summary><strong>Ground truth coverage</strong></summary>

### Ground truth coverage

<img src="images/benchmark/equal_gt_coverage_stacked_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_gt_coverage_stacked_circfl.png" width="100%"/>
<img src="images/benchmark/equal_gt_coverage_stacked_isocirc.png" width="100%"/>

circFL-seq protocol has led to the highest ground truth recovery at every tier (e.g. `discovery` R 0.50 vs CIRI-long's 0.41 and isoCirc's 0.23, from the table above); isoCirc protocol has the lowest at every tier, meaning it's reads are hardest to recover the calls from.

</details>

<details>
<summary><strong>circRNA type classification</strong></summary>

### circRNA type classification

<img src="images/benchmark/equal_type_classification_ciri_long_exon.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_ciri_long_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_ciri_long_strict_bsj.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_circfl_exon.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_circfl_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_circfl_strict_bsj.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_isocirc_exon.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_isocirc_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_isocirc_strict_bsj.png" width="100%"/>

eciRNA is the easiest type to recover on every protocol. circFL-seq protocol leads to the best recovery of every type; isoCirc protocol to the worst recovery, most severely affecting antisense circRNA (8.3% vs CIRI-long's 32.2% and circFL-seq's 36.7%).

</details>

<details>
<summary><strong>Isoform detection</strong></summary>

### Isoform detection

We are also interested to see if the callers have found the same circRNA on isoform level, which is the main point of using long-read sequencing. We evaluated their isoform detection performance with three tiers, each one requiring a BSJ boundary and an exon overlap check at different stringency levels: `full` (exact BSJ, 100% overlap, fully identical structure), `strict` (BSJ within tolerance window, 99% exon overlap, allows small boundary noise but can be easily reconstructed by follow up analyses), `relaxed` (BSJ within tolerance, 90% overlap, allowing a part of mismatched sequence).

<img src="images/benchmark/pr_scatter_isoform_ciri_long.png" width="100%"/>
<img src="images/benchmark/pr_scatter_isoform_circfl.png" width="100%"/>
<img src="images/benchmark/pr_scatter_isoform_isocirc.png" width="100%"/>

Among ground-truth loci where the same BSJ has more than one isoform (isoform family), we checked how often each caller recovers the exact isoform with different levels of stringency.

<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_isocirc.png" width="100%"/>

<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_isocirc.png" width="100%"/>

% of ground-truth isoforms called at 99% reciprocal overlap:

| Tier / tool | CIRI-long | circFL-seq | isoCirc |
| --- | --- | --- | --- |
| CIRI-long | 43.1% | 46.1% | 2.9% |
| isoCirc | 24.5% | 31.4% | 13.7% |
| circNICK-lrs | 26.5% | 27.5% | 28.4% |
| circFL-seq | 44.1% | 52.9% | 32.4% |
| `discovery` | 58.8% | 64.7% | 41.2% |
| `balanced` | 52.9% | 58.8% | 29.4% |
| `high_confidence` | 38.2% | 41.2% | 23.5% |

circFL-seq protocol allows resolution of the most isoforms at every tier. isoCirc protocol is the hardest to resolve isoform families on for every tier.

<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_isocirc.png" width="100%"/>

<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_isocirc.png" width="100%"/>

Per family plot. Shows % of isoform families fully resolved, at 99% reciprocal overlap:

| Tier / tool | CIRI-long | circFL-seq | isoCirc |
| --- | --- | --- | --- |
| CIRI-long alone | 9.8% | 17.6% | 0.0% |
| isoCirc alone | 0.0% | 5.9% | 0.0% |
| circNICK-lrs alone | 0.0% | 0.0% | 0.0% |
| circFL-seq alone | 11.8% | 25.5% | 0.0% |
| `discovery` | 29.4% | 35.3% | 11.8% |
| `balanced` | 17.6% | 23.5% | 0.0% |
| `high_confidence` | 0.0% | 0.0% | 0.0% |

Getting one isoform from the family right while missing or misclassifying others is common; fully resolving an isoform family is rare for every tier and protocol. circFL-seq protocol gives the best odds of a full family resolution at every tier; isoCirc protocol remains the hardest for all tools. Only `discovery` fully resolves any isoform families on isoCirc protocol, and circNICK-lrs and `high_confidence` never fully resolve any isoform family on any protocol.

</details>

<details>
<summary><strong>Expression quantification and reproducibility</strong></summary>

### Locus level expression quantification

<img src="images/benchmark/equal_expression_correlation_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc.png" width="100%"/>

<img src="images/benchmark/equal_expression_correlation_ciri_long_strict_bsj.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl_strict_bsj.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_strict_bsj.png" width="100%"/>

<img src="images/benchmark/equal_expression_correlation_ciri_long_exon_based.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl_exon_based.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_exon_based.png" width="100%"/>

Tool may report a correct circRNA, but misquantify it. It is important to measure, as for many research tasks that rely on differential expression precise circRNA quantification on locus and isoform level is required. We measured it with Pearson r, 
Spearman rho, R², read-count correlation against ground truth and cross-method agreement, per tool/tier. 

Pearson r by protocol (Locus level, relaxed-BSJ matching):

| Caller | CIRI-long | circFL-seq | isoCirc |
| --- | --- | --- | --- |
| CIRI-long | 0.911 | 0.935 | 0.308 |
| isoCirc | 0.528 | 0.576 | 0.403 |
| circNICK-lrs | 0.739 | 0.767 | 0.836 |
| circFL-seq | 0.765 | 0.785 | 0.694 |
| `high_confidence` | 0.948 | 0.933 | 0.802 |
| `balanced` | 0.952 | 0.936 | 0.822 |
| `discovery` | 0.911 | 0.895 | 0.814 |

isoCirc protocol is the hardest to quantify on for every caller except circNICK-lrs, whose isoCirc-protocol correlation (r=0.836) is its best of the 3 protocols, while quantification and structural-detection quality (weakest for circNICK-lrs on isoCirc protocol) are in the opposite directions. CIRI-long and circFL-seq protocol are close for every other caller (e.g. `discovery` 0.911/0.895), CIRI-long protocol usually slightly ahead. Every nanocirc tier correlates with ground truth expression more strongly than any single tool, on every protocol. isoCirc alone has the weakest quantification of any tool/tier on CIRI-long and circFL-seq protocol (r 0.53-0.58), despite it's high precision on circRNA detection. CIRI-long is the strongest quantifier among individual tools on its own and circFL-seq protocol, but the weakest on isoCirc protocol (r=0.308): the handful of calls it does make there barely correlate with the truth, worse even than isoCirc's own 0.403.

#### Isoform-level expression quantification

The table above reports relaxed-BSJ matching, which collapses several real isoforms sharing one BSJ into a single point. The tables below instead use `full`, `strict`, `relaxed` evaluation tiers, same as for [Isoform detection](#isoform-detection). 

Pearson r (n):

**`full` (exact BSJ, 100% exon overlap):**

| Caller | CIRI-long r (n) | circFL-seq r (n) | isoCirc r (n) |
| --- | --- | --- | --- |
| CIRI-long | 0.908 (7,660) | 0.931 (8,931) | 0.326 (603) |
| isoCirc | 0.530 (4,030) | 0.576 (5,314) | 0.403 (2,254) |
| circNICK-lrs | 0.615 (1,606) | 0.699 (1,838) | 0.807 (1,892) |
| circFL-seq | 0.766 (6,187) | 0.787 (7,342) | 0.705 (4,966) |
| `high_confidence` | 0.949 (7,451) | 0.953 (9,106) | 0.810 (2,652) |
| `balanced` | 0.946 (9,870) | 0.947 (12,035) | 0.831 (5,859) |
| `discovery` | 0.946 (10,124) | 0.947 (12,268) | 0.822 (6,357) |

**`strict` (BSJ within tolerance, 99% exon overlap):**

| Caller | CIRI-long r (n) | circFL-seq r (n) | isoCirc r (n) |
| --- | --- | --- | --- |
| CIRI-long | 0.909 (8,705) | 0.932 (10,701) | 0.314 (671) |
| isoCirc | 0.528 (4,053) | 0.574 (5,339) | 0.401 (2,264) |
| circNICK-lrs | 0.615 (1,696) | 0.690 (1,932) | 0.799 (1,978) |
| circFL-seq | 0.767 (6,276) | 0.785 (7,476) | 0.705 (5,060) |
| `high_confidence` | 0.949 (7,565) | 0.953 (9,274) | 0.807 (2,664) |
| `balanced` | 0.947 (10,558) | 0.948 (13,124) | 0.832 (5,927) |
| `discovery` | 0.947 (10,815) | 0.948 (13,377) | 0.822 (6,428) |

**`relaxed` (BSJ within tolerance, 90% exon overlap):**

| Caller | CIRI-long r (n) | circFL-seq r (n) | isoCirc r (n) |
| --- | --- | --- | --- |
| CIRI-long | 0.910 (9,700) | 0.933 (11,751) | 0.317 (748) |
| isoCirc | 0.528 (4,075) | 0.574 (5,379) | 0.402 (2,280) |
| circNICK-lrs | 0.618 (1,900) | 0.690 (2,141) | 0.804 (2,176) |
| circFL-seq | 0.767 (6,329) | 0.786 (7,528) | 0.692 (5,311) |
| `high_confidence` | 0.948 (7,721) | 0.951 (9,450) | 0.806 (2,691) |
| `balanced` | 0.948 (11,434) | 0.947 (13,983) | 0.831 (6,161) |
| `discovery` | 0.948 (11,706) | 0.948 (14,250) | 0.821 (6,703) |

Every nanocirc tier still correlates far more strongly than any raw tool at every tier of stringency. isoCirc protocol remains the hardest to quantify for every caller except circNICK-lrs, which quantifies best on isoCirc protocol at every tier (r 0.81 vs 0.61-0.70 on the other two); CIRI-long protocol and circFL-seq protocol stay close for everyone else.

<img src="images/benchmark/equal_expression_correlation_ciri_long_full.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl_full.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_full.png" width="100%"/>

<img src="images/benchmark/equal_expression_correlation_ciri_long_strict.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl_strict.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_strict.png" width="100%"/>

<img src="images/benchmark/equal_expression_correlation_ciri_long_relaxed.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl_relaxed.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_relaxed.png" width="100%"/>

Detection sensitivity by ground-truth expression tertile (low/mid/high read count):

| Tier | Protocol | Low | Mid | High |
| --- | --- | --- | --- | --- |
| `discovery` | CIRI-long | 34.5% | 41.0% | 50.6% |
| `discovery` | circFL-seq | 42.5% | 49.7% | 59.8% |
| `discovery` | isoCirc | 15.7% | 24.2% | 36.8% |
| `balanced` | CIRI-long | 33.2% | 40.1% | 49.6% |
| `balanced` | circFL-seq | 41.4% | 48.5% | 58.9% |
| `balanced` | isoCirc | 13.8% | 21.6% | 34.1% |
| `high_confidence` | CIRI-long | 17.1% | 22.8% | 31.6% |
| `high_confidence` | circFL-seq | 21.5% | 27.6% | 38.1% |
| `high_confidence` | isoCirc | 2.8% | 5.8% | 13.0% |

Sensitivity grows with expression at every tier and protocol (more reads means more chances to be called). `balanced` stays close to `discovery` (within 1-2 points at every tertile), while `high_confidence` sits meaningfully lower but still detects a substantial share of loci, even at low expression. circFL-seq protocol allows the highest sensitivity at every tertile and tier; isoCirc protocol the lowest, most severely at `high_confidence`, where even highly-expressed circRNAs are detected less than half as often as on circFL-seq (13.0% vs 38.1%).

<img src="images/benchmark/equal_sensitivity_by_expression_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_circfl.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_isocirc.png" width="100%"/>

Same tertiles, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison:

<img src="images/benchmark/equal_sensitivity_by_expression_ciri_long_isoform.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_circfl_isoform.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_isocirc_isoform.png" width="100%"/>

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl.png" width="100%"/>

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long_strict_bsj.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl_strict_bsj.png" width="100%"/>

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long_exon_based.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl_exon_based.png" width="100%"/>

Pairwise read-count agreement across every tool/tier combination, ground truth included. nanocirc's own tiers agree with ground truth better, than for any other caller and agree with each other almost perfectly, since they share the same remap-based quantification regardless of which tool found the locus. Agreement of raw tool calls with ground truth is looser and protocol-dependent; isoCirc shows the widest scatter against every other caller and ground truth on every protocol, consistent with its weak expression correlation above, CIRI-long shows the best agreement with ground truth among single tools on 2 out of 3 protocols.

Same reproducibility matrix, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison, so agreement is evaluated on isoform level, each real isoform scored against its own matched caller isoform (see [Isoform-level expression quantification](#isoform-level-expression-quantification) above).

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long_full.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl_full.png" width="100%"/>

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long_strict.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl_strict.png" width="100%"/>

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long_relaxed.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl_relaxed.png" width="100%"/>

### Caller overlap

Ground-truth loci recovered by each of the four tools are pooled across runs and drawn as an UpSet-style plot, where a locus counts as recovered if it was matched in at least one run. Two versions were used in the original benchmark paper: exon-based, which requires correct internal structure, and relaxed BSJ, which requires only genomic position. This is now extended to all six matching strategies used elsewhere in this benchmark: strict BSJ and exon-based (exact/structural), relaxed BSJ (coordinate-tolerant), and the full/strict/relaxed isoform tiers (structure plus BSJ tolerance, tightest to loosest).

The point of this plot is agreement, not any one tool's own recall. Of every circRNA any tool finds, only ~9% are found by all four tools together (9.2% on the species benchmark, 9.1% on the wet-lab CIRI-long protocol), and ~37-42% are found by just by one tool. Each tool is seeing circRNAs through its own lens, with no single reliable "best caller" to choose, while a naive union of all four tools' calls would inherit every tool's false positives, and a naive intersection would keep only the rare loci every tool happens to agree on, losing most of the recall. That disagreement is the core motivation for nanocirc pipeline and it's merging logic. See [Why we chose vote-based merging and not simple union/intersection approaches](methods.md#why-we-chose-vote-based-merging-and-not-simple-unionintersection-approaches) for how NanoCirc's tiers reconcile this instead.

<img src="images/benchmark/upset_plot_equal_ciri_long_strict_bsj.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl_strict_bsj.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_strict_bsj.png" width="100%"/>

<img src="images/benchmark/upset_plot_equal_ciri_long.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc.png" width="100%"/>

<img src="images/benchmark/upset_plot_equal_ciri_long_exon_based.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl_exon_based.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_exon_based.png" width="100%"/>

<img src="images/benchmark/upset_plot_equal_ciri_long_full.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl_full.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_full.png" width="100%"/>

<img src="images/benchmark/upset_plot_equal_ciri_long_strict.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl_strict.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_strict.png" width="100%"/>

<img src="images/benchmark/upset_plot_equal_ciri_long_relaxed.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_circfl_relaxed.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_relaxed.png" width="100%"/>

</details>

**Protocol results:**

- **circFL-seq protocol generates reads that work well with all tools**, at every tier: highest recall, highest F1, and the highest BSJ-confirmation rate of the 3 protocols (59.2%).
- **CIRI-long protocol reads are also compatible with all tools**, with a lower BSJ-confirmation rate (51.7%) and lower recall at every tier than circFL-seq.
- **isoCirc protocol reads are the hardest for the tools to work with**, despite carrying the exact same ~9,800 underlying circRNAs. Its recall is the lowest of the 3 protocols at every tier (BSJ-confirmation rate 41.9%), but precision remains comparable or higher to other wet-lab approaches. CIRI-long is almost incompatible with this wet lab protocol, and that was verified by running CIRI-long on real isoCIRC dataset as well - the tool never managed to finish the analysis on a full dataset.

Because the circRNA catalog, its abundances, and the simulated read amount were identical across all 3 comparisons, this difference comes entirely from how each protocol's wet-lab induced error profile and read length characteristics change the reads NanoSim generates from them.

**Recommendations**: if your protocol choice is flexible, circFL-seq or CIRI-long library prep will give nanocirc's tiers meaningfully more usable recall than isoCirc protocol. If you already generated your reads and are working with isoCirc data, you can prioritise `high_confidence` tier which reaches high precision on this read type, and expect potential incompatibility with the CIRI-long tool.

<details>
<summary><strong>Statistical significance</strong></summary>

### Statistical significance

We tested whether the protocol-to-protocol differences above are real or could be noise. Each pair matches every caller/circRNA-type/run point between the two protocols, using these statistical methods (used the same way for every "statistical significance" comparison in this document):

- **Linear regression (slope, intercept, R²)**: fits one side's value from the other's (value = intercept + slope × other value). Slope near 1 with a small intercept means the two sides are close to equal plus a fixed shift. R² is how well that line predicts the value, whatever the size of the shift.
- **Kendall's tau-b**: ranks the 8 callers best to worst on each side, then checks how well the two rankings agree.
- **Mantel-Haenszel odds ratio + Breslow-Day test**: uses the raw TP/FP/FN counts directly, stratified by circRNA type. The odds ratio is how many times more likely a call is to be correct, or a ground-truth circRNA is recovered, on one side versus the other. Breslow-Day then checks if that odds ratio is the same size for every type; p < 0.05 means the effect size depends on circRNA type.
- **Wilcoxon signed-rank**: tests if one side's values are systematically higher than the other's, without assuming a normal distribution.
- **Cohen's d (paired)**: how big that gap is, in standard deviations. ~0.2 is small, ~0.5 medium, ~0.8+ large.

#### CIRI-long vs circFL-seq

##### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | CIRI-long | circFL-seq | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 99 | 0.705 | 0.713 | +0.008 | 0.988 | 0.016 | 0.986 | 0.059 | 0.255 |
| Recall | 105 | 0.241 | 0.296 | +0.055 | 1.129 | 0.024 | 0.973 | 2×10⁻¹⁷ | 1.146 |
| F1 | 105 | 0.305 | 0.359 | +0.054 | 1.042 | 0.041 | 0.969 | 2×10⁻¹⁷ | 1.343 |

Precision fits an almost perfect 1:1 line (R² = 0.99, slope ≈ 1) and the gap itself is small and not significant (p = 0.059, d = 0.26). Recall and F1 are both significantly and substantially higher on circFL-seq (p ≈ 2×10⁻¹⁷, d ≈ 1.2-1.3), with a similarly tight fit (R² ≥ 0.97).

##### Per type, by metric (n=21 paired points per row pooling all 7 callers × 3 runs, except ciRNA precision n=18 and intergenic precision n=18: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | CIRI-long | circFL-seq | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 21 | 0.629 | 0.626 | −0.002 | 0.967 | 0.019 | 0.996 | 0.562 | −0.173 |
| eciRNA | Recall | 21 | 0.515 | 0.567 | +0.052 | 0.992 | 0.057 | 0.968 | 1×10⁻⁶ | 2.087 |
| eciRNA | F1 | 21 | 0.553 | 0.584 | +0.031 | 1.017 | 0.022 | 0.980 | 1×10⁻⁶ | 1.408 |
| EIciRNA | Precision | 21 | 0.650 | 0.670 | +0.020 | 0.995 | 0.023 | 0.965 | 0.137 | 0.407 |
| EIciRNA | Recall | 21 | 0.113 | 0.180 | +0.066 | 1.602 | −0.002 | 0.984 | 1×10⁻⁶ | 1.935 |
| EIciRNA | F1 | 21 | 0.189 | 0.277 | +0.089 | 1.445 | 0.005 | 0.984 | 1×10⁻⁶ | 2.139 |
| ciRNA | Precision | 18 | 0.882 | 0.904 | +0.022 | 0.727 | 0.263 | 0.950 | 0.006 | 0.808 |
| ciRNA | Recall | 21 | 0.118 | 0.166 | +0.048 | 1.437 | −0.003 | 0.994 | 2×10⁻⁴ | 1.480 |
| ciRNA | F1 | 21 | 0.198 | 0.264 | +0.066 | 1.325 | 0.001 | 0.996 | 2×10⁻⁴ | 1.745 |
| antisense | Precision | 21 | 0.648 | 0.639 | −0.010 | 0.982 | 0.002 | 0.999 | 0.058 | −0.584 |
| antisense | Recall | 21 | 0.157 | 0.180 | +0.023 | 1.134 | 0.002 | 0.991 | 2×10⁻⁴ | 1.041 |
| antisense | F1 | 21 | 0.245 | 0.272 | +0.027 | 1.089 | 0.005 | 0.992 | 2×10⁻⁴ | 1.075 |
| intergenic | Precision | 18 | 0.746 | 0.759 | +0.013 | 0.953 | 0.048 | 0.979 | 0.031 | 0.509 |
| intergenic | Recall | 21 | 0.300 | 0.387 | +0.087 | 1.312 | −0.006 | 0.999 | 2×10⁻⁴ | 1.083 |
| intergenic | F1 | 21 | 0.341 | 0.398 | +0.057 | 1.150 | 0.006 | 0.999 | 2×10⁻⁴ | 1.542 |

Every type's recall and F1 fit almost perfectly (R² ≥ 0.97) and circFL-seq is significantly ahead on both, for every type (p ≤ 2×10⁻⁴). Precision is smaller in magnitude and mostly not significant: only ciRNA and intergenic significantly favor circFL-seq (p ≤ 0.031); eciRNA, EIciRNA, and antisense show no significant precision difference.

##### Does caller ranking transfer?

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 7 | 0.905 (p = 0.003) | 0.964 (p = 5×10⁻⁴) |
| Recall | 7 | 1.000 (p = 4×10⁻⁴) | 1.000 (p = <10⁻³⁰⁰) |
| F1 | 7 | 1.000 (p = 4×10⁻⁴) | 1.000 (p = <10⁻³⁰⁰) |

Ranking transfers essentially perfectly. By F1:

| Caller | F1 CIRI-long | F1 circFL-seq | Rank CIRI-long | Rank circFL-seq |
| --- | --- | --- | --- | --- |
| balanced | 0.445 | 0.519 | 1 | 1 |
| discovery | 0.428 | 0.495 | 2 | 2 |
| CIRI-long | 0.375 | 0.448 | 3 | 3 |
| high_confidence | 0.329 | 0.387 | 4 | 4 |
| circFL-seq | 0.302 | 0.352 | 5 | 5 |
| isoCirc | 0.206 | 0.256 | 6 | 6 |
| circNICK-lrs | 0.050 | 0.053 | 7 | 7 |

All seven callers have the exact same F1 rank on both protocols.

##### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=9)

| Tool | Metric | n | CIRI-long | circFL-seq | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.695 | 0.733 | +0.038 | 0.913 | 0.098 | 0.925 | 1×10⁻⁴ | 1.226 |
| CIRI-long | Recall | 15 | 0.293 | 0.361 | +0.067 | 0.954 | 0.081 | 0.947 | 6×10⁻⁵ | 1.568 |
| CIRI-long | F1 | 15 | 0.375 | 0.448 | +0.073 | 0.800 | 0.147 | 0.807 | 6×10⁻⁵ | 1.329 |
| isoCirc | Precision | 15 | 0.886 | 0.873 | −0.012 | 1.110 | −0.110 | 0.929 | 0.140 | −0.371 |
| isoCirc | Recall | 15 | 0.135 | 0.178 | +0.043 | 1.218 | 0.013 | 0.997 | 6×10⁻⁵ | 1.279 |
| isoCirc | F1 | 15 | 0.206 | 0.256 | +0.051 | 1.102 | 0.030 | 0.993 | 6×10⁻⁵ | 1.855 |
| circNICK-lrs | Precision | 9 | 0.091 | 0.087 | −0.004 | 1.006 | −0.004 | 0.989 | 0.249 | −0.463 |
| circNICK-lrs | Recall | 15 | 0.057 | 0.064 | +0.008 | 1.137 | 0.000 | 1.000 | 0.028 | 0.554 |
| circNICK-lrs | F1 | 15 | 0.050 | 0.053 | +0.003 | 1.061 | 0.000 | 1.000 | 0.028 | 0.597 |
| circFL-seq | Precision | 15 | 0.679 | 0.669 | −0.009 | 0.980 | 0.004 | 0.997 | 0.095 | −0.501 |
| circFL-seq | Recall | 15 | 0.212 | 0.253 | +0.041 | 0.973 | 0.046 | 0.972 | 6×10⁻⁵ | 1.312 |
| circFL-seq | F1 | 15 | 0.302 | 0.352 | +0.050 | 0.981 | 0.056 | 0.962 | 6×10⁻⁵ | 1.196 |
| high_confidence | Precision | 15 | 0.900 | 0.895 | −0.006 | 1.123 | −0.116 | 0.980 | 0.252 | −0.384 |
| high_confidence | Recall | 15 | 0.231 | 0.282 | +0.051 | 1.058 | 0.037 | 0.989 | 6×10⁻⁵ | 2.208 |
| high_confidence | F1 | 15 | 0.329 | 0.387 | +0.058 | 0.969 | 0.068 | 0.979 | 6×10⁻⁵ | 2.047 |
| balanced | Precision | 15 | 0.741 | 0.764 | +0.023 | 0.897 | 0.100 | 0.949 | 0.008 | 0.913 |
| balanced | Recall | 15 | 0.374 | 0.463 | +0.089 | 1.149 | 0.033 | 0.963 | 6×10⁻⁵ | 1.467 |
| balanced | F1 | 15 | 0.445 | 0.519 | +0.074 | 0.932 | 0.105 | 0.973 | 6×10⁻⁵ | 2.424 |
| discovery | Precision | 15 | 0.696 | 0.718 | +0.022 | 0.938 | 0.065 | 0.955 | 0.041 | 0.657 |
| discovery | Recall | 15 | 0.381 | 0.470 | +0.089 | 1.144 | 0.034 | 0.967 | 6×10⁻⁵ | 1.482 |
| discovery | F1 | 15 | 0.428 | 0.495 | +0.067 | 0.858 | 0.128 | 0.986 | 6×10⁻⁵ | 2.345 |

Every caller's recall and F1 is significantly higher on circFL-seq (p ≤ 6×10⁻⁵) with a good fit (R² ≥ 0.81). For precision, CIRI-long, `balanced`, and `discovery` show a small but significant circFL-seq precision edge (p ≤ 0.041), while isoCirc, circNICK-lrs, circFL-seq, and `high_confidence` show no significant precision difference at all between the two protocols.

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 1.15 | 1.09–1.21 | 3×10⁻⁸ | 7×10⁻⁶ |
| CIRI-long | Recall | 5 | 1.43 | 1.38–1.49 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 0.83 | 0.74–0.92 | 4×10⁻⁴ | 0.015 |
| isoCirc | Recall | 5 | 1.50 | 1.42–1.57 | <10⁻³⁰⁰ | 7×10⁻⁴ |
| circNICK-lrs | Precision | 5 | 1.00 | 0.93–1.07 | 0.968 | 0.120 |
| circNICK-lrs | Recall | 5 | 1.19 | 1.10–1.28 | 7×10⁻⁶ | 0.801 |
| circFL-seq | Precision | 5 | 0.92 | 0.85–1.00 | 0.051 | 0.874 |
| circFL-seq | Recall | 5 | 1.32 | 1.27–1.38 | <10⁻³⁰⁰ | 3×10⁻¹⁵ |
| high_confidence | Precision | 5 | 0.90 | 0.83–0.97 | 0.009 | 0.093 |
| high_confidence | Recall | 5 | 1.38 | 1.32–1.43 | <10⁻³⁰⁰ | 6×10⁻⁵ |
| balanced | Precision | 5 | 1.10 | 1.05–1.16 | 4×10⁻⁵ | 0.002 |
| balanced | Recall | 5 | 1.60 | 1.54–1.66 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 1.07 | 1.02–1.11 | 0.003 | 6×10⁻⁵ |
| discovery | Recall | 5 | 1.61 | 1.55–1.67 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

Every caller has significantly higher odds of recovering a ground-truth circRNA (recall) on circFL-seq (OR 1.19–1.61). Precision odds are mixed: isoCirc and `high_confidence` actually favor CIRI-long (OR 0.83/0.90, p ≤ 0.009), circNICK-lrs and circFL-seq show no significant precision difference, and CIRI-long/`balanced`/`discovery` favor circFL-seq modestly (OR 1.07–1.15). 4 of 14 caller/outcome combinations have a protocol effect that's uniform across circRNA types (circNICK-lrs's precision and recall, circFL-seq's precision, `high_confidence`'s precision); the other 10 combinations' boost size still depends on which circRNA type it is.

**Conclusion.** Correlation is excellent for both metrics (aggregate R² ≥ 0.97, essentially every type and caller individually significant on recall/F1). Caller ranking transfers perfectly (Kendall's tau-b = 1.00 on recall/F1, all 7 callers same rank by F1). Precision differences are real but small and caller/type-dependent, with several callers and types showing no significant precision difference between protocols while circFL-seq is better for recall/F1.


#### CIRI-long vs isoCirc

##### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | CIRI-long | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 99 | 0.705 | 0.645 | -0.059 | 0.905 | 0.008 | 0.777 | 0.001 | -0.445 |
| Recall | 105 | 0.241 | 0.120 | -0.121 | 0.510 | -0.003 | 0.616 | 8×10⁻¹⁵ | -0.902 |
| F1 | 105 | 0.305 | 0.164 | -0.141 | 0.603 | -0.020 | 0.582 | 1×10⁻¹⁶ | -1.037 |

Precision fits reasonably well (R² = 0.78) with a moderate, significant gap (d = -0.45). Recall and F1 are both significantly and substantially lower on isoCirc protocol (p < 10⁻¹⁴, d ≈ -0.9 to -1.0), and the fit is noisier than the CIRI-long-vs-circFL-seq comparison above (R² ≈ 0.58-0.62 vs ≥ 0.97 there): isoCirc protocol doesn't just shift every caller down by a similar amount, it hits some callers, especially CIRI-long itself, far harder than others.

##### Per type, by metric (n=21 paired points per row pooling all 7 callers × 3 runs, except ciRNA precision n=18 and intergenic precision n=18: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | CIRI-long | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 21 | 0.629 | 0.642 | +0.014 | 0.900 | 0.077 | 0.787 | 0.203 | 0.136 |
| eciRNA | Recall | 21 | 0.515 | 0.315 | -0.200 | 0.098 | 0.265 | 0.011 | 1×10⁻⁵ | -1.114 |
| eciRNA | F1 | 21 | 0.553 | 0.392 | -0.160 | 0.554 | 0.086 | 0.273 | 1×10⁻⁵ | -1.057 |
| EIciRNA | Precision | 21 | 0.650 | 0.570 | -0.081 | 0.933 | -0.037 | 0.657 | 0.103 | -0.465 |
| EIciRNA | Recall | 21 | 0.113 | 0.075 | -0.038 | 0.468 | 0.022 | 0.173 | 0.024 | -0.619 |
| EIciRNA | F1 | 21 | 0.189 | 0.117 | -0.072 | 0.499 | 0.023 | 0.231 | 5×10⁻⁵ | -0.802 |
| ciRNA | Precision | 18 | 0.882 | 0.855 | -0.028 | 0.482 | 0.429 | 0.083 | 0.369 | -0.193 |
| ciRNA | Recall | 21 | 0.118 | 0.056 | -0.061 | 0.289 | 0.022 | 0.142 | 2×10⁻⁴ | -0.847 |
| ciRNA | F1 | 21 | 0.198 | 0.099 | -0.099 | 0.376 | 0.025 | 0.204 | 2×10⁻⁴ | -0.907 |
| antisense | Precision | 21 | 0.648 | 0.547 | -0.102 | 0.849 | -0.004 | 0.950 | 6×10⁻⁴ | -0.990 |
| antisense | Recall | 21 | 0.157 | 0.037 | -0.120 | 0.269 | -0.005 | 0.752 | 4×10⁻⁴ | -1.247 |
| antisense | F1 | 21 | 0.245 | 0.068 | -0.178 | 0.327 | -0.013 | 0.740 | 3×10⁻⁴ | -1.332 |
| intergenic | Precision | 18 | 0.746 | 0.644 | -0.102 | 1.112 | -0.185 | 0.801 | 0.002 | -1.043 |
| intergenic | Recall | 21 | 0.300 | 0.115 | -0.185 | 0.485 | -0.030 | 0.803 | 2×10⁻⁴ | -1.267 |
| intergenic | F1 | 21 | 0.341 | 0.143 | -0.198 | 0.509 | -0.030 | 0.581 | 2×10⁻⁴ | -1.282 |

Every type's recall and F1 is significantly lower on isoCirc protocol (p ≤ 2×10⁻⁴), same direction as the aggregate. Precision is smaller and type-dependent: antisense and intergenic are significantly lower on isoCirc protocol (p ≤ 0.002), while eciRNA, EIciRNA, and ciRNA show no significant precision difference (p ≥ 0.10).

##### Does caller ranking transfer?

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 7 | 0.810 (p = 0.011) | 0.929 (p = 0.003) |
| Recall | 7 | 0.429 (p = 0.239) | 0.536 (p = 0.215) |
| F1 | 7 | 0.429 (p = 0.239) | 0.536 (p = 0.215) |

Precision ranking transfers well (tau-b = 0.81, p = 0.011), but recall/F1 ranking does not transfer reliably (tau-b = 0.43, not significant). By F1:

| Caller | F1 CIRI-long | F1 isoCirc | Rank CIRI-long | Rank isoCirc |
| --- | --- | --- | --- | --- |
| balanced | 0.445 | 0.302 | 1 | 1 |
| discovery | 0.428 | 0.284 | 2 | 2 |
| CIRI-long | 0.375 | 0.042 | 3 | 7 |
| high_confidence | 0.329 | 0.107 | 4 | 5 |
| circFL-seq | 0.302 | 0.235 | 5 | 3 |
| isoCirc | 0.206 | 0.124 | 6 | 4 |
| circNICK-lrs | 0.050 | 0.052 | 7 | 6 |

The two nanocirc tiers that lead on CIRI-long protocol (`balanced`, `discovery`) still lead on isoCirc protocol. Everything else reshuffles: the CIRI-long tool itself falls from 3rd to 7th, since it barely functions on isoCirc-protocol reads (see below), while circFL-seq and isoCirc both move up a few places because they degrade less than the callers around them.

##### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=9)

| Tool | Metric | n | CIRI-long | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.695 | 0.661 | -0.033 | 0.368 | 0.405 | 0.157 | 0.561 | -0.278 |
| CIRI-long | Recall | 15 | 0.293 | 0.022 | -0.271 | 0.132 | -0.016 | 0.909 | 6×10⁻⁵ | -1.673 |
| CIRI-long | F1 | 15 | 0.375 | 0.042 | -0.333 | 0.338 | -0.085 | 0.844 | 6×10⁻⁵ | -3.962 |
| isoCirc | Precision | 15 | 0.886 | 0.889 | +0.003 | 0.934 | 0.062 | 0.807 | 0.778 | 0.067 |
| isoCirc | Recall | 15 | 0.135 | 0.076 | -0.060 | 0.630 | -0.010 | 0.967 | 6×10⁻⁵ | -1.049 |
| isoCirc | F1 | 15 | 0.206 | 0.124 | -0.081 | 0.697 | -0.019 | 0.951 | 6×10⁻⁵ | -1.197 |
| circNICK-lrs | Precision | 9 | 0.091 | 0.077 | -0.014 | 0.949 | -0.009 | 0.952 | 0.028 | -0.820 |
| circNICK-lrs | Recall | 15 | 0.057 | 0.066 | +0.009 | 1.148 | 0.001 | 0.999 | 0.028 | 0.605 |
| circNICK-lrs | F1 | 15 | 0.050 | 0.052 | +0.002 | 1.033 | 0.000 | 0.999 | 0.046 | 0.530 |
| circFL-seq | Precision | 15 | 0.679 | 0.447 | -0.231 | 0.677 | -0.012 | 0.808 | 1×10⁻⁴ | -1.441 |
| circFL-seq | Recall | 15 | 0.212 | 0.171 | -0.042 | 0.709 | 0.020 | 0.955 | 0.060 | -0.680 |
| circFL-seq | F1 | 15 | 0.302 | 0.235 | -0.067 | 0.732 | 0.014 | 0.977 | 2×10⁻⁴ | -1.057 |
| high_confidence | Precision | 15 | 0.900 | 0.836 | -0.065 | 0.967 | -0.035 | 0.245 | 0.421 | -0.512 |
| high_confidence | Recall | 15 | 0.231 | 0.069 | -0.163 | 0.549 | -0.058 | 0.778 | 6×10⁻⁵ | -1.682 |
| high_confidence | F1 | 15 | 0.329 | 0.107 | -0.223 | 0.685 | -0.119 | 0.676 | 6×10⁻⁵ | -2.005 |
| balanced | Precision | 15 | 0.741 | 0.735 | -0.005 | 0.798 | 0.144 | 0.311 | 0.762 | -0.041 |
| balanced | Recall | 15 | 0.374 | 0.213 | -0.161 | 0.521 | 0.019 | 0.684 | 6×10⁻⁵ | -1.194 |
| balanced | F1 | 15 | 0.445 | 0.302 | -0.143 | 0.601 | 0.035 | 0.458 | 6×10⁻⁵ | -1.038 |
| discovery | Precision | 15 | 0.696 | 0.645 | -0.051 | 1.115 | -0.132 | 0.819 | 0.026 | -0.618 |
| discovery | Recall | 15 | 0.381 | 0.221 | -0.160 | 0.593 | -0.005 | 0.787 | 6×10⁻⁵ | -1.335 |
| discovery | F1 | 15 | 0.428 | 0.284 | -0.144 | 0.689 | -0.011 | 0.562 | 6×10⁻⁵ | -1.288 |

`high_confidence` has full n=15 coverage on precision. CIRI-long's own recall collapses from 29% on its own protocol to 2.2% on isoCirc protocol (d = -1.67), and its F1 drops by an order of magnitude (0.375 → 0.042, d = -3.96, the single largest effect in this table): it is barely functioning on isoCirc-protocol reads. Every other caller's recall and F1 is also significantly lower on isoCirc protocol (p ≤ 0.060), consistent with the aggregate result, but nowhere near as severely. Precision is smaller and mostly not significant: CIRI-long, isoCirc, `high_confidence`, and `balanced` show no significant precision difference, while circNICK-lrs, circFL-seq, and `discovery` do.

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 1.07 | 0.93–1.22 | 0.345 | 7×10⁻⁹ |
| CIRI-long | Recall | 5 | 0.04 | 0.04–0.05 | <10⁻³⁰⁰ | 3×10⁻¹⁶ |
| isoCirc | Precision | 5 | 0.95 | 0.83–1.08 | 0.422 | 0.618 |
| isoCirc | Recall | 5 | 0.47 | 0.44–0.50 | <10⁻³⁰⁰ | 6×10⁻⁶ |
| circNICK-lrs | Precision | 5 | 0.91 | 0.85–0.98 | 0.013 | 1×10⁻⁴ |
| circNICK-lrs | Recall | 5 | 1.23 | 1.14–1.32 | 9×10⁻⁸ | 0.462 |
| circFL-seq | Precision | 5 | 0.27 | 0.25–0.29 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circFL-seq | Recall | 5 | 0.73 | 0.70–0.76 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| high_confidence | Precision | 5 | 1.08 | 0.96–1.22 | 0.226 | 4×10⁻¹⁰ |
| high_confidence | Recall | 5 | 0.18 | 0.17–0.20 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced | Precision | 5 | 0.98 | 0.92–1.03 | 0.388 | <10⁻³⁰⁰ |
| balanced | Recall | 5 | 0.41 | 0.39–0.42 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 0.78 | 0.74–0.82 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Recall | 5 | 0.40 | 0.39–0.42 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

CIRI-long's recall odds ratio (0.04, 95% CI 0.04–0.05) is the largest single effect measured anywhere in this document: its odds of recovering a ground-truth circRNA are roughly 25× lower on isoCirc-protocol reads than on its own protocol's reads. Every caller's recall odds are significantly lower on isoCirc protocol (OR 0.04–0.73), except circNICK-lrs, which is actually higher (OR = 1.23): it is the one tool this protocol change doesn't hurt. Precision odds are smaller and mixed: circNICK-lrs, circFL-seq, and `discovery` favor CIRI-long protocol, while CIRI-long's own precision, isoCirc's, `high_confidence`'s, and `balanced`'s show no significant precision difference. Only 1 of 14 caller/outcome combinations has a protocol effect that's uniform across circRNA types (isoCirc's precision), the least homogeneity of any comparison in this document; every other combination's boost size depends heavily on which circRNA type it is.

**Conclusion.** Correlation is reasonable for precision (R² = 0.78) but weaker for recall/F1 (R² ≈ 0.58–0.62), and unlike the CIRI-long-vs-circFL-seq comparison above, caller ranking does not transfer reliably on recall/F1 (Kendall's tau-b = 0.43, not significant), though it does on precision (tau-b = 0.81). The reason is CIRI-long itself: its recall and F1 on isoCirc-protocol reads are close to zero (recall OR = 0.04, F1 Cohen's d = -3.96), dragging it from 3rd to 7th place in the ranking and confirming numerically what the recommendations at the top of this page already state qualitatively, that CIRI-long is close to incompatible with isoCirc-protocol reads. Every other caller is also significantly weaker on recall/F1 on isoCirc protocol, but far less severely, and circNICK-lrs is actually unaffected or slightly better (recall OR = 1.23).

#### circFL-seq vs isoCirc

##### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | circFL-seq | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 99 | 0.713 | 0.645 | -0.067 | 0.917 | -0.008 | 0.790 | 1×10⁻⁴ | -0.521 |
| Recall | 105 | 0.296 | 0.120 | -0.176 | 0.433 | -0.009 | 0.583 | 7×10⁻¹⁷ | -1.079 |
| F1 | 105 | 0.359 | 0.164 | -0.195 | 0.549 | -0.033 | 0.541 | 2×10⁻¹⁷ | -1.292 |

Same pattern as the CIRI-long comparison above, larger in size: circFL-seq protocol scores the highest of the 3 protocols to begin with, so the drop to isoCirc-protocol numbers is the biggest gap measured in this document for recall/F1 (d ≈ -1.1 to -1.3). Precision again fits well (R² = 0.79) with a moderate, significant gap (d = -0.52).

##### Per type, by metric (n=21 paired points per row pooling all 7 callers × 3 runs, except ciRNA precision n=18 and intergenic precision n=18: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | circFL-seq | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 21 | 0.626 | 0.642 | +0.016 | 0.941 | 0.053 | 0.807 | 0.128 | 0.171 |
| eciRNA | Recall | 21 | 0.567 | 0.315 | -0.252 | 0.069 | 0.276 | 0.005 | 1×10⁻⁵ | -1.376 |
| eciRNA | F1 | 21 | 0.584 | 0.392 | -0.191 | 0.541 | 0.077 | 0.275 | 1×10⁻⁶ | -1.250 |
| EIciRNA | Precision | 21 | 0.670 | 0.570 | -0.100 | 0.928 | -0.052 | 0.668 | 0.026 | -0.588 |
| EIciRNA | Recall | 21 | 0.180 | 0.075 | -0.105 | 0.248 | 0.030 | 0.127 | 1×10⁻⁵ | -1.210 |
| EIciRNA | F1 | 21 | 0.277 | 0.117 | -0.160 | 0.310 | 0.031 | 0.189 | 5×10⁻⁶ | -1.358 |
| ciRNA | Precision | 18 | 0.904 | 0.855 | -0.050 | 0.935 | 0.009 | 0.173 | 0.766 | -0.382 |
| ciRNA | Recall | 21 | 0.166 | 0.056 | -0.109 | 0.182 | 0.026 | 0.116 | 2×10⁻⁴ | -1.098 |
| ciRNA | F1 | 21 | 0.264 | 0.099 | -0.165 | 0.272 | 0.027 | 0.188 | 2×10⁻⁴ | -1.197 |
| antisense | Precision | 21 | 0.639 | 0.547 | -0.092 | 0.858 | -0.002 | 0.936 | 0.002 | -0.847 |
| antisense | Recall | 21 | 0.180 | 0.037 | -0.143 | 0.246 | -0.007 | 0.820 | 3×10⁻⁴ | -1.281 |
| antisense | F1 | 21 | 0.272 | 0.068 | -0.204 | 0.310 | -0.017 | 0.792 | 2×10⁻⁴ | -1.386 |
| intergenic | Precision | 18 | 0.759 | 0.644 | -0.115 | 1.207 | -0.272 | 0.875 | 1×10⁻³ | -1.383 |
| intergenic | Recall | 21 | 0.387 | 0.115 | -0.272 | 0.371 | -0.029 | 0.810 | 2×10⁻⁴ | -1.235 |
| intergenic | F1 | 21 | 0.398 | 0.143 | -0.255 | 0.435 | -0.030 | 0.562 | 2×10⁻⁴ | -1.372 |

Every type's recall and F1 is significantly lower on isoCirc protocol (p ≤ 2×10⁻⁴), same direction as the aggregate and as the CIRI-long comparison above. Precision is smaller and type-dependent: EIciRNA, antisense, and intergenic are significantly lower on isoCirc protocol (p ≤ 0.026), while eciRNA and ciRNA show no significant precision difference.

##### Does caller ranking transfer?

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 7 | 0.905 (p = 0.003) | 0.964 (p = 5×10⁻⁴) |
| Recall | 7 | 0.429 (p = 0.239) | 0.536 (p = 0.215) |
| F1 | 7 | 0.429 (p = 0.239) | 0.536 (p = 0.215) |

Precision ranking transfers strongly (tau-b = 0.91, p = 0.003); recall and F1 ranking do not (p = 0.239), the same pattern as the CIRI-long comparison above. By F1:

| Caller | F1 circFL-seq | F1 isoCirc | Rank circFL-seq | Rank isoCirc |
| --- | --- | --- | --- | --- |
| balanced | 0.519 | 0.302 | 1 | 1 |
| discovery | 0.495 | 0.284 | 2 | 2 |
| CIRI-long | 0.448 | 0.042 | 3 | 7 |
| high_confidence | 0.387 | 0.107 | 4 | 5 |
| circFL-seq | 0.352 | 0.235 | 5 | 3 |
| isoCirc | 0.256 | 0.124 | 6 | 4 |
| circNICK-lrs | 0.053 | 0.052 | 7 | 6 |

Same reshuffle as the CIRI-long comparison above, and for the same reason: CIRI-long falls from 3rd to 7th because it barely functions on isoCirc-protocol reads regardless of which protocol it's being compared against, while circFL-seq and isoCirc both move up a few places.

##### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=9)

| Tool | Metric | n | circFL-seq | isoCirc | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.733 | 0.661 | -0.071 | 0.383 | 0.381 | 0.153 | 0.041 | -0.608 |
| CIRI-long | Recall | 15 | 0.361 | 0.022 | -0.338 | 0.131 | -0.025 | 0.867 | 6×10⁻⁵ | -2.129 |
| CIRI-long | F1 | 15 | 0.448 | 0.042 | -0.406 | 0.346 | -0.113 | 0.701 | 6×10⁻⁵ | -5.304 |
| isoCirc | Precision | 15 | 0.873 | 0.889 | +0.016 | 0.852 | 0.145 | 0.893 | 0.182 | 0.398 |
| isoCirc | Recall | 15 | 0.178 | 0.076 | -0.103 | 0.516 | -0.016 | 0.965 | 6×10⁻⁵ | -1.160 |
| isoCirc | F1 | 15 | 0.256 | 0.124 | -0.132 | 0.626 | -0.036 | 0.936 | 6×10⁻⁵ | -1.473 |
| circNICK-lrs | Precision | 9 | 0.087 | 0.077 | -0.010 | 0.954 | -0.006 | 0.985 | 0.028 | -1.039 |
| circNICK-lrs | Recall | 15 | 0.064 | 0.066 | +0.002 | 1.009 | 0.001 | 0.999 | 0.046 | 0.502 |
| circNICK-lrs | F1 | 15 | 0.053 | 0.052 | -0.001 | 0.974 | 0.000 | 1.000 | 0.116 | -0.449 |
| circFL-seq | Precision | 15 | 0.669 | 0.447 | -0.222 | 0.690 | -0.015 | 0.811 | 6×10⁻⁵ | -1.427 |
| circFL-seq | Recall | 15 | 0.253 | 0.171 | -0.082 | 0.720 | -0.012 | 0.961 | 1×10⁻³ | -1.426 |
| circFL-seq | F1 | 15 | 0.352 | 0.235 | -0.117 | 0.727 | -0.021 | 0.964 | 6×10⁻⁵ | -1.762 |
| high_confidence | Precision | 15 | 0.895 | 0.836 | -0.059 | 0.827 | 0.096 | 0.230 | 0.421 | -0.459 |
| high_confidence | Recall | 15 | 0.282 | 0.069 | -0.213 | 0.493 | -0.070 | 0.711 | 6×10⁻⁵ | -1.870 |
| high_confidence | F1 | 15 | 0.387 | 0.107 | -0.280 | 0.653 | -0.146 | 0.590 | 6×10⁻⁵ | -2.276 |
| balanced | Precision | 15 | 0.764 | 0.735 | -0.029 | 0.925 | 0.029 | 0.354 | 0.524 | -0.226 |
| balanced | Recall | 15 | 0.463 | 0.213 | -0.250 | 0.415 | 0.021 | 0.597 | 6×10⁻⁵ | -1.396 |
| balanced | F1 | 15 | 0.519 | 0.302 | -0.217 | 0.623 | -0.022 | 0.440 | 6×10⁻⁵ | -1.603 |
| discovery | Precision | 15 | 0.718 | 0.645 | -0.073 | 1.149 | -0.180 | 0.801 | 0.026 | -0.831 |
| discovery | Recall | 15 | 0.470 | 0.221 | -0.249 | 0.487 | -0.008 | 0.717 | 6×10⁻⁵ | -1.525 |
| discovery | F1 | 15 | 0.495 | 0.284 | -0.211 | 0.811 | -0.118 | 0.582 | 6×10⁻⁵ | -2.097 |

`high_confidence` again has full n=15 coverage on precision. CIRI-long's F1 drop (0.448 → 0.042, d = -5.30) is the single largest effect in this table, and its recall alone drops from 36% to 2.2% (d = -2.13). Every other caller's recall and F1 is also significantly lower on isoCirc protocol (p ≤ 1×10⁻³), same direction as the CIRI-long-vs-isoCirc comparison above but generally with larger gaps, since circFL-seq protocol is the strongest baseline of the 3. `circNICK-lrs`'s F1 is the only metric with no significant difference at all (p = 0.116), consistent with it being the one caller this protocol change barely touches. Precision is again smaller and mixed: isoCirc, `high_confidence`, and `balanced` show no significant precision-side effect, while CIRI-long, circFL-seq, and `discovery` do.

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 0.96 | 0.84–1.10 | 0.556 | 2×10⁻¹⁰ |
| CIRI-long | Recall | 5 | 0.03 | 0.03–0.03 | <10⁻³⁰⁰ | 2×10⁻¹⁶ |
| isoCirc | Precision | 5 | 1.13 | 1.00–1.27 | 0.060 | 0.039 |
| isoCirc | Recall | 5 | 0.31 | 0.30–0.33 | <10⁻³⁰⁰ | 2×10⁻⁴ |
| circNICK-lrs | Precision | 5 | 0.92 | 0.85–0.98 | 0.012 | 0.019 |
| circNICK-lrs | Recall | 5 | 1.03 | 0.96–1.11 | 0.402 | 0.303 |
| circFL-seq | Precision | 5 | 0.28 | 0.26–0.30 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circFL-seq | Recall | 5 | 0.56 | 0.54–0.59 | <10⁻³⁰⁰ | 4×10⁻¹² |
| high_confidence | Precision | 5 | 1.22 | 1.08–1.37 | 1×10⁻³ | 9×10⁻¹¹ |
| high_confidence | Recall | 5 | 0.13 | 0.13–0.14 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced | Precision | 5 | 0.86 | 0.82–0.91 | 9×10⁻⁸ | <10⁻³⁰⁰ |
| balanced | Recall | 5 | 0.26 | 0.25–0.27 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 0.72 | 0.68–0.75 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Recall | 5 | 0.26 | 0.25–0.27 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

CIRI-long's recall odds ratio (0.03, 95% CI 0.03–0.03) is again the largest effect in this table, in line with the CIRI-long-vs-isoCirc comparison above (OR = 0.04 there): its odds of recovering a ground-truth circRNA are roughly 30× lower on isoCirc-protocol reads than on circFL-seq protocol's reads. Every caller's recall odds are significantly lower on isoCirc protocol (OR 0.03–0.56) except circNICK-lrs, which shows no significant difference (OR = 1.03, p = 0.402): the same tool the CIRI-long comparison also found unaffected. Precision odds are mixed: circFL-seq, circNICK-lrs, `balanced`, and `discovery` favor circFL-seq protocol, `high_confidence` favors isoCirc protocol, and isoCirc's own precision shows no significant difference. Only 1 of 14 caller/outcome combinations has a protocol effect that's uniform across circRNA types (circNICK-lrs's recall), consistent with the CIRI-long comparison's low homogeneity count.

**Conclusion.** Same overall picture as the CIRI-long-vs-isoCirc comparison above, scaled up: circFL-seq protocol is the strongest baseline of the 3, so its gap to isoCirc-protocol numbers is the largest measured in this document for recall/F1 (d ≈ -1.1 to -1.3 aggregate, up to d = -5.30 for CIRI-long's own F1). Precision ranking transfers strongly (tau-b = 0.91) but recall/F1 ranking does not (tau-b = 0.43, not significant), again because CIRI-long collapses on isoCirc-protocol reads (recall OR = 0.03) regardless of which other protocol it's being compared against. circNICK-lrs is consistently the one caller isoCirc protocol doesn't hurt.

</details>

<details>
<summary><strong>What are these false positives?</strong></summary>

### What are these false positives?

We traced every FP-supporting read back to the simulator's own ground truth (single-tool discovery exon-based FPs, all 3 runs per protocol, mean % ± SD across runs), separately for each protocol's own reads:

#### CIRI-long protocol

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.1 ± 0.0 |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 86.8 ± 0.4 | 76.8 ± 0.7 | 93.7 ± 0.5 | 2.5 ± 0.3 |
| Real circ read, correct strand, pure boundary noise | 1.2 ± 0.2 | 4.9 ± 1.0 | 0.0 ± 0.0 | 31.1 ± 0.7 |
| Real circ read, but circRNA reported at a different/distant locus | 12.0 ± 0.4 | 18.3 ± 0.7 | 6.3 ± 0.5 | 66.3 ± 0.5 |
| FP locus overlaps a genomic repeat | 18.3 ± 0.2 | 31.2 ± 1.0 | 19.9 ± 0.4 | 62.5 ± 0.3 |
| n FP loci (pooled, 3 runs) | 3,566 | 1,660 | 513 | 7,999 |
| n supporting reads (pooled, 3 runs) | 60,498 | 14,942 | 5,274 | 27,927 |

#### circFL-seq protocol

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.1 ± 0.0 |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 84.3 ± 0.9 | 73.1 ± 2.1 | 93.1 ± 0.5 | 2.4 ± 0.4 |
| Real circ read, correct strand, pure boundary noise | 2.7 ± 0.6 | 7.1 ± 2.5 | 0.0 ± 0.0 | 23.7 ± 1.3 |
| Real circ read, but circRNA reported at a different/distant locus | 13.0 ± 0.6 | 19.8 ± 0.4 | 6.9 ± 0.4 | 73.8 ± 1.7 |
| FP locus overlaps a genomic repeat | 20.6 ± 0.8 | 37.8 ± 0.8 | 28.7 ± 1.5 | 58.7 ± 0.2 |
| n FP loci (pooled, 3 runs) | 3,702 | 2,136 | 635 | 9,316 |
| n supporting reads (pooled, 3 runs) | 48,061 | 20,115 | 7,227 | 39,033 |

#### isoCirc protocol

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.2 ± 0.1 |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 72.1 ± 8.8 | 28.9 ± 0.2 | 91.8 ± 0.9 | 1.9 ± 0.1 |
| Real circ read, correct strand, pure boundary noise | 0.7 ± 0.6 | 60.2 ± 0.3 | 0.0 ± 0.0 | 18.5 ± 1.2 |
| Real circ read, but circRNA reported at a different/distant locus | 27.3 ± 8.5 | 10.9 ± 0.3 | 8.2 ± 0.9 | 79.4 ± 1.3 |
| FP locus overlaps a genomic repeat | 11.2 ± 1.6 | 52.4 ± 0.8 | 22.4 ± 1.4 | 57.4 ± 0.5 |
| n FP loci (pooled, 3 runs) | 241 | 2,269 | 354 | 10,995 |
| n supporting reads (pooled, 3 runs) | 592 | 13,832 | 1,103 | 29,882 |

Many FP come from `antisense`-type circRNA, where the simulator deliberately sets the circRNA's strand opposite to its host gene; since strand agreement is required by every matching strategy, no boundary-tolerance setting can rescue these. That pattern is present for CIRI-long and isoCirc on all 3 protocols (72-94% wrong-strand), and for circFL-seq on the CIRI-long and circFL-seq protocols (73-77%), but circFL-seq's FP composition changes on isoCirc-protocol reads: wrong-strand drops to 28.9% and "correct strand, pure boundary noise" goes up to 60.2%, the opposite of its own pattern everywhere else and of every other tool on this protocol. CIRI-long's own FP composition on isoCirc protocol is also unusual: barely any single-tool FP loci at all (241, an order of magnitude below its 3,566-3,702 on the other 2 protocols, consistent with how few calls it makes here overall), and its second-largest bucket is "different/distant locus" (27.3%) rather than boundary noise, which it almost never has (0.7%). circNICK-lrs's FPs are structurally different from the other 3 tools throughout: mostly boundary noise and distant-locus calls (66-79%), rarely a strand flip (2-3%), and consistently the most repeat-associated (57-63% vs 11-38% for the other tools).

**Practical implication**: most `discovery` false positives still are at the right gene/circRNA, just with imprecise boundary or wrong strand. That makes them usable for gene- or locus-level expression signal analysis, but we recommend to re-verify this tier's low confidence calls in wet-lab, for example by Sanger sequencing. Do not blindly trust their reported BSJ/structure/strand and deposit them into the databases.

</details>

<details>
<summary><strong>What discovery still misses?</strong></summary>

### What discovery still misses?

We also characterized every ground-truth circRNA against whether any tool recovered it at all, per protocol (missed + recovered = 100% of the ~9,861-circRNA GT catalog, pooled across 3 runs):

#### CIRI-long protocol

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 60.6% | 39.4% |
| Mean mature length (bp) | 3,156 | 353 |
| Mean rolling-circle copy number | 7.8 | 8.7 |
| Median simulated read count | 49 | 71 |
| Canonical GT/AG splice site at BSJ | 63.2% | 76.7% |
| No gene overlap at all (intergenic) | 21.8% | 14.6% |
| Falls in a gene with ≥5 isoforms | 72.1% | 78.2% |
| Overlaps a genomic repeat | 49.7% | 54.0% |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 19.8% | 80.2% |
| `antisense` | 66.3% | 33.7% |
| `ciRNA` (purely intronic) | 68.6% | 31.4% |
| `EIciRNA` (exon+intron) | 79.3% | 20.7% |
| `intergenic` | 69.4% | 30.6% |

#### circFL-seq protocol

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 52.0% | 48.0% |
| Mean mature length (bp) | 3,502 | 477 |
| Mean rolling-circle copy number | 7.5 | 8.9 |
| Median simulated read count | 48 | 69 |
| Canonical GT/AG splice site at BSJ | 62.0% | 75.6% |
| No gene overlap at all (intergenic) | 22.6% | 14.9% |
| Falls in a gene with ≥5 isoforms | 71.2% | 78.1% |
| Overlaps a genomic repeat | 49.3% | 53.7% |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 13.4% | 86.7% |
| `antisense` | 61.4% | 38.6% |
| `ciRNA` (purely intronic) | 56.6% | 43.4% |
| `EIciRNA` (exon+intron) | 67.0% | 33.0% |
| `intergenic` | 62.5% | 37.5% |

#### isoCirc protocol

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 76.3% | 23.7% |
| Mean mature length (bp) | 2,509 | 572 |
| Mean rolling-circle copy number | 7.6 | 9.9 |
| Median simulated read count | 49 | 97 |
| Canonical GT/AG splice site at BSJ | 63.6% | 84.6% |
| No gene overlap at all (intergenic) | 21.6% | 10.3% |
| Falls in a gene with ≥5 isoforms | 72.2% | 82.0% |
| Overlaps a genomic repeat | 49.7% | 57.1% |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 41.0% | 59.0% |
| `antisense` | 89.9% | 10.1% |
| `ciRNA` (purely intronic) | 84.8% | 15.2% |
| `EIciRNA` (exon+intron) | 78.8% | 21.2% |
| `intergenic` | 87.8% | 12.2% |

circFL-seq protocol has the largest recovery of the 3 (48.0% of the GT catalog found by at least one caller), CIRI-long protocol is close behind (39.4%), and isoCirc protocol has by far the smallest (23.7%, barely more than half of CIRI-long protocol's). The missed-vs-recovered pattern is consistent across all 3 protocols regardless of that gap: missed circRNAs are longer (~4-9x), lower-expressed, and use non-canonical splice sites more often than recovered ones. `eciRNA` is always the easiest biotype to recover in all wet-lab protocols. On isoCirc protocol every biotype's recall collapses further than on the other 2 protocols, most severely for `antisense` (10.1% recovered, versus 33.7% on CIRI-long protocol and 38.6% on circFL-seq protocol), also due to the same wrong-strand FP pattern documented above rather than due to a different reason. Repeat overlap is not a meaningful driver of recovery on any protocol.

</details>

<details>
<summary><strong>isoCirc protocol: with CIRI-long off by default</strong></summary>

#### isoCirc protocol: with CIRI-long off by default

`--wet_lab isocirc` defaults `--run_cirilong` to `false` (see [Which tools work with my wet-lab protocol and species?](#which-tools-work-with-my-wet-lab-protocol-and-species) above). The isoCirc protocol column above used all 4 tools to estimate real performance of each tool; but we show here the results on the same dataset with CIRI-long left out, next to the all-4-tools numbers for comparison. 

<img src="images/benchmark/pr_scatter_equal_isocirc_no_cirilong.png" width="100%"/>

Mean precision (P) / recall (R) / F1 across the same 3 runs, exon-based matching, plus each individual tool's own raw calls:

| Tier / tool | CIRI-long on | CIRI-long off (new default) |
| --- | --- | --- |
| isoCirc alone | P 0.84 / R 0.08 / F1 0.14 | P 0.84 / R 0.08 / F1 0.14 |
| circFL-seq alone | P 0.58 / R 0.19 / F1 0.29 | P 0.58 / R 0.19 / F1 0.29 |
| CIRI-long alone | P 0.77 / R 0.03 / F1 0.05 | not run |
| circNICK-lrs alone | P 0.20 / R 0.10 / F1 0.13 | P 0.20 / R 0.10 / F1 0.13 |
| `discovery` | P 0.59 / R 0.24 / F1 0.34 | P 0.59 / R 0.24 / F1 0.34 |
| `balanced` | P 0.73 / R 0.23 / F1 0.35 | P 0.73 / R 0.23 / F1 0.35 |
| `high_confidence` | P 0.82 / R 0.07 / F1 0.13 | P 0.82 / R 0.07 / F1 0.13 |

CIRI-long was the weakest tool on this protocol (F1 0.05, far behind its performance on the CIRI-long/circFL-seq protocols).

Dropping CIRI-long costs essentially nothing at any tier in this benchmark: P/R/F1 are identical with CIRI-long on vs off. Against the real-world risk documented above (CIRI-long hanging for 24+ days on real isoCirc-protocol data), this is why CIRI-long off is now the default for this protocol. If you wish to try to run it on isoCirc data despite the risk, you can turn it on yourself with `--run_cirilong true`.

circNICK-lrs wet-lab protocol reads were not tested, as all 3 other tools expect rolling-circle generated data with multiple copies of circRNA per read, making circNICK-lrs protocol compatible only with circNICK-lrs tool, thus it wouldn't benefit from nanocirc's consensus voting algorithm and it defeats the purpose of the pipeline. You still can run circNICK-lrs inside the nanocirc pipeline and benefit from nanocirc's upstream patches to circNICK-lrs tool and nanocirc quantification approach, but without confidence scoring or a sensible Precision-Recall curve.

The plots below cover this default configuration (CIRI-long excluded) in the same depth as the other two protocols get in the main section above, since this is the tool set an actual `--wet_lab isocirc` run uses. The main section's own isoCirc-protocol plots instead use a with-CIRI-long tool set, for direct comparability with the CIRI-long protocol and circFL-seq protocol columns.

<img src="images/benchmark/equal_gt_coverage_stacked_isocirc_no_cirilong.png" width="100%"/>

<img src="images/benchmark/equal_type_classification_isocirc_exon_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_isocirc_relaxed_bsj_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_type_classification_isocirc_strict_bsj_no_cirilong.png" width="100%"/>

<img src="images/benchmark/pr_scatter_isoform_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_isocirc_no_cirilong.png" width="100%"/>

<img src="images/benchmark/equal_expression_correlation_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_strict_bsj_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_exon_based_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_full_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_strict_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc_relaxed_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_isocirc_isoform_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_strict_bsj_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_exon_based_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_full_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_strict_no_cirilong.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc_relaxed_no_cirilong.png" width="100%"/>

### Caller agreement

Ground-truth loci recovered by each of the four tools are pooled across runs and drawn as an UpSet-style plot, where a locus counts as recovered if it was matched in at least one run. Two versions were used in the original benchmark paper: exon-based, which requires correct internal structure, and relaxed BSJ, which requires only genomic position. This is now extended to all six matching strategies used elsewhere in this benchmark: strict BSJ and exon-based (exact/structural), relaxed BSJ (coordinate-tolerant), and the full/strict/relaxed isoform tiers (structure plus BSJ tolerance, tightest to loosest).

The point of this plot is agreement, not any one tool's own recall. Of every circRNA any tool finds, only about 9% are found by all four tools together (9.2% on the species benchmark, 9.1% on the wet-lab CIRI-long protocol), and roughly 40% are found by exactly one tool. Each tool is effectively viewing circRNAs through its own lens, with no single reliable "ground truth caller" to defer to. That disagreement is the core motivation for NanoCirc: a naive union of all four tools' calls would inherit every tool's false positives, and a naive intersection would keep only the rare loci every tool happens to agree on, trading away most of the recall. See [Why we chose vote-based merging and not simple union/intersection approaches](methods.md#why-we-chose-vote-based-merging-and-not-simple-unionintersection-approaches) for how NanoCirc's tiers reconcile this instead.

<img src="images/benchmark/upset_plot_equal_isocirc_strict_bsj_no_cirilong.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_no_cirilong.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_exon_based_no_cirilong.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_full_no_cirilong.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_strict_no_cirilong.png" width="100%"/>
<img src="images/benchmark/upset_plot_equal_isocirc_relaxed_no_cirilong.png" width="100%"/>

</details>

## Species benchmark: human vs mouse

Wet-lab protocol was fixed to CIRI-long for both species; each species used its own genome, GTF, circRNA database and in-silico generated by database fit circRNA catalog. Any difference below comes from species circRNA biology and genome/annotation/circRNA database differences, not the wet-lab protocol.

<details>
<summary><strong>Structural precision/recall</strong></summary>

### Structural precision/recall

<img src="images/benchmark/pr_scatter_human.png" width="100%"/>
<img src="images/benchmark/pr_scatter_mouse.png" width="100%"/>

Mean precision (P) / recall (R) / F1 across 3 runs, exon-based matching, per nanocirc tier:

| Tier | Human | Mouse |
| --- | --- | --- |
| `discovery` | P 0.67 / R 0.36 / F1 0.47 | P 0.67 / R 0.49 / F1 0.56 |
| `balanced` | P 0.76 / R 0.35 / F1 0.48 | P 0.79 / R 0.48 / F1 0.60 |
| `high_confidence` | P 0.86 / R 0.20 / F1 0.32 | P 0.90 / R 0.29 / F1 0.43 |

</details>

<details>
<summary><strong>Ground truth coverage</strong></summary>

### Ground truth coverage

<img src="images/benchmark/gt_coverage_stacked_human.png" width="100%"/>
<img src="images/benchmark/gt_coverage_stacked_mouse.png" width="100%"/>

Mouse has higher ground truth recovery than human at every tier (e.g. `discovery` R 0.49 vs human's 0.36, from the table above), the same species gap the summary table shows.

</details>

<details>
<summary><strong>circRNA type classification</strong></summary>

### circRNA type classification

<img src="images/benchmark/type_classification_human_exon.png" width="100%"/>
<img src="images/benchmark/type_classification_human_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_exon.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_strict_bsj.png" width="100%"/>

Same TP/FP vs GT breakdown as the GT coverage plot above, split by circRNA biotype (`eciRNA`, `EIciRNA`, `ciRNA`, `antisense`, `intergenic`, see [docs/methods.md](methods.md#circrna-type-classification)). `discovery`-tier recovery per type, isoform-level matching:

| Species | eciRNA | EIciRNA | ciRNA | antisense | intergenic |
| --- | --- | --- | --- | --- | --- |
| Human | 75.2% | 17.7% | 23.3% | 27.9% | 17.4% |
| Mouse | 82.5% | 40.8% | 44.7% | 45.3% | 15.6% |

Mouse recovers every type better than human except intergenic, where human is slightly ahead (17.4% vs 15.6%), the one type/species combination where human derived circRNAs win. eciRNA is the easiest type on both species; intergenic is hardest on both.

</details>

<details>
<summary><strong>Isoform detection</strong></summary>

### Isoform detection

The structural plots show if a specific part of circRNA structure was found at its locus. But we are also interested to see if the callers have found the exactly same molecule, especially on isoform level, which is the main point of using long-read sequencing. We evaluated that with three tiers, each pairing a BSJ boundary check with an exon overlap check: `full` (exact BSJ, 100% overlap, byte identical structure), `strict` (BSJ within tolerance window, 99% exon overlap, allows little boundary noise), `relaxed` (BSJ within tolerance, 90% overlap, allowing a part of mismatched sequence).

<img src="images/benchmark/pr_scatter_isoform_human.png" width="100%"/>
<img src="images/benchmark/pr_scatter_isoform_mouse.png" width="100%"/>

Long-read protocols are used mainly to resolve internal exon structure, not just BSJ position. Among ground-truth loci where the same BSJ has more than one real isoform (isoform family), we checked how often each caller recovers the exact isoform, versus calls the right BSJ with the wrong structure, versus misses the BSJ entirely.

<img src="images/benchmark/isoform_detection_reciprocal_overlap_99_human.png" width="100%"/>
<img src="images/benchmark/isoform_detection_reciprocal_overlap_99_mouse.png" width="100%"/>

<img src="images/benchmark/isoform_detection_reciprocal_overlap_90_human.png" width="100%"/>
<img src="images/benchmark/isoform_detection_reciprocal_overlap_90_mouse.png" width="100%"/>

Green+grey sum to 100% of the isoform-family ground truth; red shows the calls that found the right BSJ but didn't match internal structure criteria. 99% reciprocal exon overlap allows only trivial boundary noise; 90% allows a small sequence mismatch while still requiring the same structure, not a different one. % of ground-truth isoforms called exactly right, at 99% reciprocal overlap:

| Tier / tool | Human | Mouse |
| --- | --- | --- |
| CIRI-long alone | 18.2% | 30.8% |
| isoCirc alone | 12.8% | 20.8% |
| circNICK-lrs alone | 3.1% | 5.6% |
| circFL-seq alone | 17.8% | 23.2% |
| `discovery` | 29.6% | 41.9% |
| `balanced` | 28.0% | 40.6% |
| `high_confidence` | 14.7% | 21.4% |

Mouse resolves more isoforms than human at every tier and tool, consistent with mouse's higher overall recall in this benchmark. Both species are far below the wet-lab-protocol benchmark's exact-isoform rates above (e.g. `discovery` 29.6%/41.9% here vs 58.8-64.7% on CIRI-long/circFL-seq protocol), since this dataset uses a more isoform-family-rich catalog. CIRI-long and circFL-seq lead the raw tools on both species (18-31%); circNICK-lrs is the weakest by a wide margin (3.1%/5.6%), and `high_confidence` clears every raw tool on both species too.

<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_mouse.png" width="100%"/>

<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_90_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_90_mouse.png" width="100%"/>

Per family plot instead of per isoform: did the caller get every isoform in the family right, some of them, or none. % of isoform families fully resolved, at 99% reciprocal overlap:

| Tier / tool | Human | Mouse |
| --- | --- | --- |
| CIRI-long alone | 5.7% | 13.7% |
| isoCirc alone | 1.9% | 10.3% |
| circNICK-lrs alone | 0.0% | 0.0% |
| circFL-seq alone | 8.6% | 15.0% |
| `discovery` | 13.4% | 26.0% |
| `balanced` | 11.5% | 23.5% |
| `high_confidence` | 4.3% | 5.4% |

Getting one isoform right while missing or misclassifying others is the most common outcome; fully resolving an isoform family is rare for every tier on both species, and `high_confidence` manages it least often (4.3%/5.4%). Mouse is consistently ahead of human here too. circNICK-lrs never fully resolves a family on either species (0.0%); circFL-seq leads the raw tools (8.6%/15.0%), ahead of `high_confidence` on human (8.6% vs 4.3%).

</details>

<details>
<summary><strong>Expression quantification and reproducibility</strong></summary>

### Expression quantification and reproducibility

<img src="images/benchmark/expression_correlation_human.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse.png" width="100%"/>

<img src="images/benchmark/expression_correlation_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_strict_bsj.png" width="100%"/>

<img src="images/benchmark/expression_correlation_human_exon_based.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_exon_based.png" width="100%"/>

Read-count correlation against ground truth (right panel) and cross-method agreement (left panel), per tool/tier. Each cell shows three metrics:

- **Pearson r**: linear correlation between predicted and true counts. Sensitive to whether counts scale proportionally, not just in the same direction.
- **Spearman rho**: rank correlation. Only checks whether higher-expressed circRNAs get higher counts, regardless of exact scale. More resistant to outliers and non-linear (but still monotonic) relationships.
- **R²** (Pearson r squared): how much of the variance in true expression is captured by the predicted counts, on a 0-1 scale.

For differential expression, Spearman rho matters to rank circRNAs correctly relative to each other; Pearson r/R² matter to ensure the fold-changes DESeq2 estimates between conditions are quantitatively trustworthy, not just directionally right. Pearson r by species (relaxed-BSJ matching):

| Caller | Human | Mouse |
| --- | --- | --- |
| CIRI-long | 0.911 | 0.906 |
| isoCirc | 0.552 | 0.667 |
| circNICK-lrs | 0.763 | 0.814 |
| circFL-seq | 0.785 | 0.805 |
| `discovery` | 0.923 | 0.881 |
| `balanced` | 0.946 | 0.899 |
| `high_confidence` | 0.955 | 0.896 |

isoCirc has the weakest quantification of any tool/tier on both species (r 0.55-0.67), same as on every wet-lab protocol above. Every nanocirc tier correlates more strongly with ground truth than any raw tool on both species (r ≥ 0.88). Unlike detection, where mouse leads, quantification quality is close between species and `high_confidence` is actually slightly stronger on human (0.955 vs 0.896).

#### Isoform-level expression quantification

The table above reports relaxed-BSJ matching, which collapses a locus with several real isoforms into a single point. The tables below instead score each real GT isoform on its own: `full`, `strict`, `relaxed` are the same 3 tiers used for [Isoform detection](#isoform-detection-1) and the isoform PR scatter above (`full` = exact BSJ + 100% exon overlap, `strict` = BSJ within tolerance + 99% overlap, `relaxed` = BSJ within tolerance + 90% overlap). A GT isoform's predicted count is whichever caller isoform structurally matched it under that tier, or 0 if none did, so n here is a genuinely stricter re-match, not the table above's n expanded, and can be smaller. Pearson r (n):

**`full` (exact BSJ, 100% exon overlap):**

| Caller | Human r (n) | Mouse r (n) |
| --- | --- | --- |
| CIRI-long | 0.908 (5,748) | 0.903 (7,522) |
| isoCirc | 0.557 (3,342) | 0.675 (6,454) |
| circNICK-lrs | 0.537 (1,515) | 0.760 (1,600) |
| circFL-seq | 0.792 (4,830) | 0.811 (6,654) |
| `discovery` | 0.952 (7,891) | 0.920 (11,260) |
| `balanced` | 0.953 (7,622) | 0.921 (11,005) |
| `high_confidence` | 0.959 (5,623) | 0.937 (9,010) |

**`strict` (BSJ within tolerance, 99% exon overlap):**

| Caller | Human r (n) | Mouse r (n) |
| --- | --- | --- |
| CIRI-long | 0.909 (6,518) | 0.904 (9,181) |
| isoCirc | 0.557 (3,363) | 0.673 (6,545) |
| circNICK-lrs | 0.538 (1,596) | 0.765 (1,685) |
| circFL-seq | 0.792 (4,892) | 0.809 (6,753) |
| `discovery` | 0.954 (8,405) | 0.922 (12,246) |
| `balanced` | 0.954 (8,130) | 0.924 (11,982) |
| `high_confidence` | 0.959 (5,699) | 0.935 (9,241) |

**`relaxed` (BSJ within tolerance, 90% exon overlap):**

| Caller | Human r (n) | Mouse r (n) |
| --- | --- | --- |
| CIRI-long | 0.909 (7,217) | 0.906 (9,824) |
| isoCirc | 0.557 (3,384) | 0.673 (6,588) |
| circNICK-lrs | 0.545 (1,732) | 0.762 (1,868) |
| circFL-seq | 0.791 (4,935) | 0.809 (6,791) |
| `discovery` | 0.953 (9,029) | 0.921 (12,561) |
| `balanced` | 0.953 (8,745) | 0.923 (12,286) |
| `high_confidence` | 0.959 (5,783) | 0.934 (9,338) |

n grows from `full` to `relaxed` on both species as the overlap requirement loosens, and mouse's n grows a bit faster throughout, consistent with mouse's richer isoform-family catalog documented in [Isoform detection](#isoform-detection-1) above. circNICK-lrs is the caller most affected by tightening the requirement: its correlation drops against the locus-level table above (r 0.54-0.76 here vs r 0.76-0.81 there) and its n is by far the smallest of any caller (1,515-1,868), since its structure calls carry real boundary noise a genuine per-isoform re-match does not forgive. Every nanocirc tier still correlates far more strongly than any raw tool at every tier of stringency, and human/mouse stay close to each other throughout, same as the table above.

<img src="images/benchmark/expression_correlation_human_full.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_full.png" width="100%"/>

<img src="images/benchmark/expression_correlation_human_strict.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_strict.png" width="100%"/>

<img src="images/benchmark/expression_correlation_human_relaxed.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_relaxed.png" width="100%"/>

<img src="images/benchmark/sensitivity_by_expression_human.png" width="100%"/>
<img src="images/benchmark/sensitivity_by_expression_mouse.png" width="100%"/>

Same tertiles, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison:

<img src="images/benchmark/sensitivity_by_expression_human_isoform.png" width="100%"/>
<img src="images/benchmark/sensitivity_by_expression_mouse_isoform.png" width="100%"/>

Detection sensitivity broken down by ground-truth expression tertile: low-expressed circRNAs are structurally harder to detect correctly than high-expressed ones, for every tool.

| Tier | Species | Low | Mid | High |
| --- | --- | --- | --- | --- |
| `discovery` | Human | 16.2% | 38.7% | 56.1% |
| `discovery` | Mouse | 26.4% | 55.3% | 68.7% |
| `balanced` | Human | 14.4% | 37.7% | 55.2% |
| `balanced` | Mouse | 24.1% | 54.7% | 68.3% |
| `high_confidence` | Human | 5.3% | 20.2% | 34.0% |
| `high_confidence` | Mouse | 11.5% | 32.1% | 44.9% |

Mouse is more sensitive than human at every tertile and tier shown, widest at `discovery`'s low tertile (26.4% vs 16.2%, a 63% relative gap) and narrowest at `high_confidence`'s high tertile, where the species gap has mostly closed.

<img src="images/benchmark/reproducibility_pairs_human.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse.png" width="100%"/>

<img src="images/benchmark/reproducibility_pairs_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_strict_bsj.png" width="100%"/>

<img src="images/benchmark/reproducibility_pairs_human_exon_based.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_exon_based.png" width="100%"/>

Pairwise read-count agreement across every tool/tier combination, ground truth included. As with the wet-lab-protocol benchmark above, nanocirc's own tiers agree with each other almost perfectly (diagonal points in every tier-vs-tier panel) since they share the same remap-based quantification; agreement against raw tool calls is looser, isoCirc again showing the widest scatter, consistent with its weaker Pearson r above.

Same reproducibility matrix, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison, so agreement is evaluated on isoform level, each real isoform scored against its own matched caller isoform (see [Isoform-level expression quantification](#isoform-level-expression-quantification-1) above).

<img src="images/benchmark/reproducibility_pairs_human_full.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_full.png" width="100%"/>

<img src="images/benchmark/reproducibility_pairs_human_strict.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_strict.png" width="100%"/>

<img src="images/benchmark/reproducibility_pairs_human_relaxed.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_relaxed.png" width="100%"/>

### Caller agreement

Ground-truth loci recovered by each of the four tools are pooled across runs and drawn as an UpSet-style plot, where a locus counts as recovered if it was matched in at least one run.

<img src="images/benchmark/upset_plot_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse_strict_bsj.png" width="100%"/>

<img src="images/benchmark/upset_plot_human.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse.png" width="100%"/>

<img src="images/benchmark/upset_plot_human_exon_based.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse_exon_based.png" width="100%"/>

<img src="images/benchmark/upset_plot_human_full.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse_full.png" width="100%"/>

<img src="images/benchmark/upset_plot_human_strict.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse_strict.png" width="100%"/>

<img src="images/benchmark/upset_plot_human_relaxed.png" width="100%"/>
<img src="images/benchmark/upset_plot_mouse_relaxed.png" width="100%"/>

</details>

**What we see:**

- **The species effect is real for recall, not for precision.** Pooled across every caller and circRNA type, mouse recall/F1 is significantly higher (Wilcoxon p < 10⁻¹⁵), but precision shows no significant species gap at all (p = 0.44), and actually reverses for antisense and intergenic circRNAs, where human precision is significantly higher (see [Statistical significance](#statistical-significance-1) below). NanoCirc's own merged tiers still score higher on mouse at all 3 tiers, on both precision and recall (table above): the F1 gap is largest at `balanced` (F1 0.48 human vs 0.60 mouse) and smallest at `discovery` (F1 0.47 vs 0.56), but that tier-level edge doesn't reflect a uniform species effect once you look past nanocirc's own merged output.

- **Quantification differences between species are not as big**, on both precision and recall. Quantification performance is close between species and `high_confidence` is actually slightly stronger on human (0.955 vs 0.896).

**Recommendations**: species-level difference for detection are present, especially for recall and F1. Published tier-level precision/recall for human or mouse are not as guaranteed to transfer across other species. Re-benchmark on your own species/catalog using nanocirc-simulator and nanocirc-benchmark and turn off circNICK-lrs for any species but human or mouse (since it's incompatible with any other species).

<details>
<summary><strong>Statistical significance</strong></summary>

### Statistical significance

Same statistical methods as the [wet-lab-protocol comparison](#statistical-significance) above (linear regression, Kendall's tau-b, Mantel-Haenszel + Breslow-Day, Wilcoxon signed-rank, Cohen's d), matching every caller/circRNA-type/run point between human and mouse.

#### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | Human | Mouse | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 96 | 0.723 | 0.742 | +0.019 | 0.942 | 0.061 | 0.843 | 0.079 | 0.195 |
| Recall | 105 | 0.210 | 0.305 | +0.096 | 1.182 | 0.057 | 0.910 | 6×10⁻¹⁷ | 1.191 |
| F1 | 105 | 0.274 | 0.378 | +0.104 | 1.119 | 0.072 | 0.909 | 2×10⁻¹⁶ | 1.381 |

Recall and F1 are well predicted by human's value (R² = 0.91) with a slope above 1: mouse's advantage is higher for callers that are stronger, not growing uniformly for all callers. Precision fits well (R² = 0.84, slope near 1) but shows no significant species gap at all (Wilcoxon p = 0.079): mouse is not reliably more precise than human once every caller and type is pooled together.

#### Per type, by metric (n=21 paired points per row pooling all 7 callers × 3 runs, except ciRNA precision n=18, antisense precision n=18, and intergenic precision n=18: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | Human | Mouse | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 21 | 0.647 | 0.754 | +0.107 | 1.076 | 0.058 | 0.913 | 1×10⁻⁶ | 1.468 |
| eciRNA | Recall | 21 | 0.476 | 0.556 | +0.081 | 1.041 | 0.061 | 0.893 | 1×10⁻⁶ | 1.695 |
| eciRNA | F1 | 21 | 0.536 | 0.635 | +0.099 | 1.209 | −0.013 | 0.964 | 1×10⁻⁶ | 2.164 |
| EIciRNA | Precision | 21 | 0.645 | 0.667 | +0.023 | 1.046 | −0.007 | 0.848 | 0.517 | 0.214 |
| EIciRNA | Recall | 21 | 0.099 | 0.199 | +0.100 | 2.169 | −0.016 | 0.918 | 1×10⁻⁶ | 1.955 |
| EIciRNA | F1 | 21 | 0.168 | 0.301 | +0.133 | 1.937 | −0.024 | 0.909 | 1×10⁻⁵ | 1.878 |
| ciRNA | Precision | 18 | 0.862 | 0.924 | +0.062 | 0.585 | 0.420 | 0.973 | 8×10⁻⁶ | 1.766 |
| ciRNA | Recall | 21 | 0.095 | 0.196 | +0.101 | 1.908 | 0.014 | 0.980 | 2×10⁻⁴ | 1.774 |
| ciRNA | F1 | 21 | 0.165 | 0.305 | +0.141 | 1.646 | 0.034 | 0.976 | 2×10⁻⁴ | 2.069 |
| antisense | Precision | 18 | 0.766 | 0.714 | −0.052 | 0.938 | −0.005 | 0.971 | 9×10⁻⁴ | −0.862 |
| antisense | Recall | 21 | 0.134 | 0.232 | +0.098 | 1.645 | 0.011 | 0.873 | 5×10⁻⁴ | 0.953 |
| antisense | F1 | 21 | 0.215 | 0.316 | +0.101 | 1.381 | 0.020 | 0.914 | 5×10⁻⁴ | 1.019 |
| intergenic | Precision | 18 | 0.719 | 0.659 | −0.060 | 0.649 | 0.192 | 0.883 | 0.007 | −0.848 |
| intergenic | Recall | 21 | 0.244 | 0.343 | +0.100 | 1.511 | −0.025 | 0.993 | 5×10⁻⁴ | 0.822 |
| intergenic | F1 | 21 | 0.287 | 0.333 | +0.047 | 1.184 | −0.006 | 0.991 | 4×10⁻⁴ | 1.002 |

Mouse is significantly ahead on recall and F1 for every type (p ≤ 9×10⁻⁴), including intergenic this time. For precision the species relationship flips by type: eciRNA and ciRNA favor mouse (p ≤ 8×10⁻⁶), EIciRNA shows no significant difference (p = 0.52), and antisense and intergenic actually favor human (negative differences, p ≤ 0.007, the only rows in this whole table where human is significantly ahead of mouse).

#### Does caller ranking transfer?

We averaged each caller's human and mouse values across type and run first, ranked the 7 callers by that average on each side, then compared the two rankings with Kendall's tau-b:

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 7 | 0.905 (p = 0.003) | 0.964 (p = 5×10⁻⁴) |
| Recall | 7 | 1.000 (p = 4×10⁻⁴) | 1.000 (p = <10⁻³⁰⁰) |
| F1 | 7 | 1.000 (p = 4×10⁻⁴) | 1.000 (p = <10⁻³⁰⁰) |

Ranking transfers significantly on every metric, and recall/F1 transfer perfectly. By F1:

| Caller | F1 human | F1 mouse | Rank human | Rank mouse |
| --- | --- | --- | --- | --- |
| balanced | 0.411 | 0.555 | 1 | 1 |
| discovery | 0.397 | 0.525 | 2 | 2 |
| CIRI-long | 0.320 | 0.442 | 3 | 3 |
| high_confidence | 0.281 | 0.400 | 4 | 4 |
| circFL-seq | 0.267 | 0.353 | 5 | 5 |
| isoCirc | 0.187 | 0.318 | 6 | 6 |
| circNICK-lrs | 0.055 | 0.055 | 7 | 7 |

All seven callers reach the exact same F1 rank on both species.

#### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=6)

F1 here is lower than the [structural precision/recall](#structural-precisionrecall-1) table above for the same tier (e.g. `discovery` 0.40/0.53 here vs 0.47/0.56 there): that table pools all circRNA types together per run before averaging across the 3 runs, while this one averages precision/recall/F1 per circRNA type first, then across type × run (15 points). The two are different aggregations of the same underlying calls, not different data.

| Tool | Metric | n | Human | Mouse | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.688 | 0.791 | +0.103 | 0.379 | 0.530 | 0.329 | 0.004 | 1.052 |
| CIRI-long | Recall | 15 | 0.241 | 0.338 | +0.097 | 1.053 | 0.084 | 0.832 | 0.002 | 1.116 |
| CIRI-long | F1 | 15 | 0.320 | 0.442 | +0.122 | 1.132 | 0.080 | 0.710 | 0.001 | 1.142 |
| isoCirc | Precision | 15 | 0.870 | 0.801 | −0.069 | 0.964 | −0.038 | 0.501 | 0.026 | −0.817 |
| isoCirc | Recall | 15 | 0.120 | 0.236 | +0.115 | 1.484 | 0.057 | 0.892 | 6×10⁻⁵ | 1.198 |
| isoCirc | F1 | 15 | 0.187 | 0.318 | +0.131 | 1.271 | 0.081 | 0.929 | 6×10⁻⁵ | 1.550 |
| circNICK-lrs | Precision | 6 | 0.166 | 0.129 | −0.037 | 4.259 | −0.578 | 0.917 | 0.438 | −0.712 |
| circNICK-lrs | Recall | 15 | 0.056 | 0.059 | +0.003 | 1.036 | 0.001 | 0.999 | 0.028 | 0.767 |
| circNICK-lrs | F1 | 15 | 0.055 | 0.055 | −0.001 | 1.036 | −0.003 | 0.993 | 0.600 | −0.096 |
| circFL-seq | Precision | 15 | 0.670 | 0.665 | −0.005 | 1.011 | −0.012 | 0.974 | 0.972 | −0.086 |
| circFL-seq | Recall | 15 | 0.184 | 0.253 | +0.069 | 1.057 | 0.059 | 0.932 | 0.002 | 1.352 |
| circFL-seq | F1 | 15 | 0.267 | 0.353 | +0.086 | 1.042 | 0.075 | 0.905 | 0.002 | 1.216 |
| high_confidence | Precision | 15 | 0.904 | 0.904 | −0.001 | 0.516 | 0.437 | 0.139 | 0.978 | −0.006 |
| high_confidence | Recall | 15 | 0.192 | 0.282 | +0.090 | 1.040 | 0.082 | 0.994 | 6×10⁻⁵ | 5.963 |
| high_confidence | F1 | 15 | 0.281 | 0.400 | +0.119 | 0.950 | 0.133 | 0.986 | 6×10⁻⁵ | 5.042 |
| balanced | Precision | 15 | 0.736 | 0.800 | +0.064 | 0.676 | 0.303 | 0.502 | 0.008 | 0.796 |
| balanced | Recall | 15 | 0.332 | 0.481 | +0.148 | 1.102 | 0.114 | 0.901 | 6×10⁻⁵ | 1.876 |
| balanced | F1 | 15 | 0.411 | 0.555 | +0.143 | 0.897 | 0.186 | 0.993 | 6×10⁻⁵ | 6.147 |
| discovery | Precision | 15 | 0.690 | 0.734 | +0.044 | 0.924 | 0.096 | 0.807 | 0.026 | 0.603 |
| discovery | Recall | 15 | 0.342 | 0.490 | +0.147 | 1.090 | 0.116 | 0.912 | 6×10⁻⁵ | 1.922 |
| discovery | F1 | 15 | 0.397 | 0.525 | +0.128 | 0.813 | 0.202 | 0.964 | 6×10⁻⁵ | 3.132 |

`high_confidence` has full n=15 coverage on precision. Every caller is significantly higher on mouse for recall and F1 (p ≤ 0.002). Precision is smaller and mixed: CIRI-long, isoCirc, `balanced`, and `discovery` show a significant species gap (p ≤ 0.026, all favoring mouse except isoCirc), while circNICK-lrs, circFL-seq, and `high_confidence` show no significant precision difference between species.

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 2.16 | 2.03–2.30 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| CIRI-long | Recall | 5 | 1.79 | 1.72–1.87 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 0.68 | 0.61–0.76 | 5×10⁻¹² | <10⁻³⁰⁰ |
| isoCirc | Recall | 5 | 2.72 | 2.59–2.86 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circNICK-lrs | Precision | 5 | 0.89 | 0.83–0.96 | 0.003 | <10⁻³⁰⁰ |
| circNICK-lrs | Recall | 5 | 1.07 | 0.99–1.15 | 0.089 | 0.335 |
| circFL-seq | Precision | 5 | 1.16 | 1.05–1.27 | 0.002 | <10⁻³⁰⁰ |
| circFL-seq | Recall | 5 | 1.66 | 1.59–1.74 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| high_confidence | Precision | 5 | 1.41 | 1.28–1.56 | 9×10⁻¹² | <10⁻³⁰⁰ |
| high_confidence | Recall | 5 | 1.81 | 1.73–1.90 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced | Precision | 5 | 1.46 | 1.38–1.54 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced | Recall | 5 | 2.08 | 2.00–2.17 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 1.22 | 1.16–1.28 | 8×10⁻¹⁶ | <10⁻³⁰⁰ |
| discovery | Recall | 5 | 2.09 | 2.00–2.17 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

Every caller has significantly higher odds of recovering a ground-truth circRNA (recall) on mouse than on human (OR 1.1–3.1), and most also have higher odds of a correct call (precision), though isoCirc's precision actually favors human (OR = 0.68, p = 5×10⁻¹²) and circNICK-lrs's precision also favors human (OR = 0.89, p = 0.003). Only 1 of the 14 caller/outcome combinations (circNICK-lrs's recall) has a species effect that's stable across circRNA types (Breslow-Day p = 0.34); every other combination's effect size depends on circRNA type.

**Conclusion.** Correlation is good for recall and F1 (R² = 0.91, perfect rank transfer), but weaker and inconsistent for precision, where several callers and types show no significant species difference at all, and antisense/intergenic precision (and isoCirc's own precision) actually favor human. Caller ranking transfers strongly (Kendall's tau-b ≥ 0.91 on every metric, perfect on recall/F1). The size of the species gap depends on circRNA type (1 of 14 Mantel-Haenszel/Breslow-Day combinations homogeneous). So the recall/F1 species gap exists, the caller ranking is essentially fully preserved between the species, and precision across species is mostly caller- and type-dependent.

</details>

<details>
<summary><strong>What are these false positives?</strong></summary>

### What are these false positives?

We traced every FP-supporting read back to the simulator's own ground truth (single-tool discovery exon-based FPs, all 3 human runs, mean % ± SD across runs):

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 87.1 ± 0.2 | 77.6 ± 2.2 | 92.4 ± 0.8 | 5.6 ± 0.5 |
| Real circ read, correct strand, pure boundary noise | 2.2 ± 0.5 | 6.2 ± 2.5 | 0.0 ± 0.0 | 34.5 ± 2.4 |
| Real circ read, but circRNA reported at a different/distant locus | 10.7 ± 0.4 | 16.3 ± 0.5 | 7.6 ± 0.8 | 59.8 ± 2.7 |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.1 ± 0.0 |
| FP locus overlaps a genomic repeat | 18.9 ± 0.1 | 34.0 ± 1.3 | 24.5 ± 3.0 | 64.6 ± 0.6 |
| n FP loci (pooled, 3 runs) | 2,679 | 1,248 | 375 | 6,844 |
| n supporting reads (pooled, 3 runs) | 63,109 | 15,341 | 4,675 | 29,099 |

Many FP come from `antisense`-type circRNA, where the simulator deliberately sets the circRNA's strand opposite to its host gene. NanoSim also mimics cDNA based RCRT-sequencing, meaning it simulates both strand sequences: one as an actual circRNA read and the other as cDNA read. The tools correctly find the genomic position of circRNA but select the host gene's (wrong) strand instead of flipping it to the antisense one (correct). Since strand agreement is required by all three matching strategies (exon-based, relaxed BSJ, strict BSJ), no boundary-tolerance setting can rescue these FP and we didn't find an easy way to override this gene-strand bias without introducing new FP elsewhere. These realistic strand misses could also happen on real data and we suggest to keep that in mind as you move onto wet-lab validation.

**Practical implication**: most `discovery` false positives still point at the right gene/circRNA locus, just with an imprecise boundary or wrong strand. That makes them usable for gene- or locus-level expression signal analysis, but not for isoform-level analysis or primer/probe design, both of which need to get the exact strand/boundary.

</details>

<details>
<summary><strong>What discovery still misses?</strong></summary>

### What discovery still misses?

We also characterized every ground-truth circRNA against whether any tool recovered it at all (mean % ± SD across 3 human runs; missed + recovered = 100% of the GT catalog, ~8,700 circRNAs/run):

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 64.8 ± 0.0 | 35.2 ± 0.0 |
| Median mature length (bp) | 1,093 ± 10 | 277 ± 2 |
| Mean rolling-circle copy number | 7.7 ± 0.1 | 8.6 ± 0.0 |
| Median simulated read count | 35 ± 0 | 99 ± 0 |
| Canonical GT/AG splice site at BSJ | 63.8 ± 0.1 | 80.2 ± 0.1 |
| No gene overlap at all (intergenic) | 19.9 ± 0.2 | 9.1 ± 0.2 |
| Falls in a gene with ≥5 isoforms | 73.7 ± 0.1 | 83.0 ± 0.5 |
| Overlaps a genomic repeat | 51.1 ± 0.5 | 57.7 ± 0.9 |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 22.9 ± 0.2 | 77.1 ± 0.2 |
| `antisense` | 68.9 ± 0.1 | 31.1 ± 0.1 |
| `ciRNA` (purely intronic) | 74.5 ± 0.2 | 25.5 ± 0.2 |
| `EIciRNA` (exon+intron) | 82.1 ± 0.1 | 17.9 ± 0.1 |
| `intergenic` | 81.3 ± 0.2 | 18.7 ± 0.2 |

Missed circRNAs are generally longer (~4x), lower-expressed (~2.8x fewer reads), and use non-canonical splice sites more often. Rolling-circle copy number is also slightly lower, meaning weaker RCRT/RCA tandem-repeat evidence for the tools that rely on it. `eciRNA` seems to be the easier type to recover; every other biotype, especially without a clean match to an annotated exon, is more likely to be missed. Repeat overlap isn't a meaningful driver of circRNA discovery.

</details>

---

## Database benchmark: annotated circRNA recovery

Other benchmarks use synthetic circRNAs generated by nanocirc-simulator and processed into reads by NanoSim. This benchmark evaluates performance using 2,000 database-derived circRNAs per type instead, sampled randomly from the human circBase/circAtlas union catalog. It uses the same NanoSim read simulation, expression distribution, detection pipeline, and quantification as every other dataset here. It shows the recovery of already-annotated circRNAs versus in silico de novo ones. It is important to note that these databases can contain circRNAs derived from short-read sequencing computational predictions that were never verified in a wet lab, are biased toward what tools can already detect (since that is how they were catalogued in the first place), and that some tools use the circRNA database itself as direct input for detection. It is expected that tools will perform better on database circRNAs, and there is no guarantee that every circRNA in the database is biologically true. Even with these limitations, it is a useful comparator and the closest approximation to biological ground truth we have.

<details>
<summary><strong>Structural precision/recall</strong></summary>

### Structural precision/recall

<img src="images/benchmark/pr_scatter_database_human.png" width="100%"/>

Mean precision (P) / recall (R) / F1 across 3 runs, exon-based matching, per nanocirc tier:

| Tier | Database (human) |
| --- | --- |
| `discovery` | P 0.79 / R 0.48 / F1 0.60 |
| `balanced` | P 0.91 / R 0.47 / F1 0.62 |
| `high_confidence` | P 0.97 / R 0.36 / F1 0.52 |

</details>

<details>
<summary><strong>Ground truth coverage</strong></summary>

### Ground truth coverage

<img src="images/benchmark/gt_coverage_stacked_database_human.png" width="100%"/>

Green (recovered) plus grey (missed) sum to 100% of ground truth; red stacks separately on top as false calls, scaled the same way.

</details>

<details>
<summary><strong>circRNA type classification</strong></summary>

### circRNA type classification

<img src="images/benchmark/type_classification_database_human_isoform_level.png" width="100%"/>

Same TP/FP vs GT breakdown as the ground truth coverage plot above, split by circRNA biotype (`eciRNA`, `EIciRNA`, `ciRNA`, `antisense`, `intergenic`, see [docs/methods.md](methods.md#circrna-type-classification)) instead of one bar per caller. `discovery`-tier recovery per type, isoform-level matching:

| Type | eciRNA | EIciRNA | ciRNA | antisense | intergenic |
| --- | --- | --- | --- | --- | --- |
| Database (human) | 79.1% | 23.9% | 48.9% | 43.0% | 33.7% |

eciRNA, ciRNA, and antisense are recovered well above the same tier's synthetic species-benchmark rates (see [circRNA type classification](#circrna-type-classification-1) above); EIciRNA and intergenic remain the two hardest types here too, the same pattern seen throughout this document.

</details>

<details>
<summary><strong>Isoform detection</strong></summary>

### Isoform detection

<img src="images/benchmark/pr_scatter_isoform_database_human.png" width="100%"/>
<img src="images/benchmark/isoform_detection_reciprocal_overlap_99_database_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_database_human.png" width="100%"/>

**Note: only 2 shared-BSJ multi-isoform families exist in this 2,000-circRNAs-per-type database sample**, versus 102 in the wet-lab-protocol benchmark and 422/912 in the species benchmark. The numbers below describe those 2 families; treat them accordingly. % of ground-truth isoforms called right, at 99% reciprocal overlap:

| Tier / tool | Isoform | Family |
| --- | --- | --- |
| CIRI-long alone | 41.7% | 33.3% |
| isoCirc alone | 33.3% | 0.0% |
| circNICK-lrs alone | 0.0% | 0.0% |
| circFL-seq alone | 50.0% | 0.0% |
| `discovery` | 75.0% | 50.0% |
| `balanced` | 75.0% | 50.0% |
| `high_confidence` | 25.0% | 0.0% |

</details>

<details>
<summary><strong>Expression quantification and reproducibility</strong></summary>

### Expression quantification and reproducibility

<img src="images/benchmark/expression_correlation_database_human.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human.png" width="100%"/>

<img src="images/benchmark/expression_correlation_database_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human_strict_bsj.png" width="100%"/>

<img src="images/benchmark/expression_correlation_database_human_exon_based.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human_exon_based.png" width="100%"/>

Read-count correlation against ground truth, and pairwise cross-tool/tier agreement, per caller. Pearson r (relaxed-BSJ matching):

| Caller | Database (human) |
| --- | --- |
| CIRI-long | 0.916 |
| isoCirc | 0.551 |
| circNICK-lrs | 0.776 |
| circFL-seq | 0.778 |
| `discovery` | 0.923 |
| `balanced` | 0.953 |
| `high_confidence` | 0.954 |

isoCirc is again the weakest quantifier, and every nanocirc tier again correlates more strongly with ground truth than any raw tool, the same pattern as the wet-lab-protocol and species benchmarks above. Detection sensitivity by ground-truth expression tertile:

| Tier | Low | Mid | High |
| --- | --- | --- | --- |
| `discovery` | 28.8% | 53.9% | 62.4% |
| `balanced` | 26.0% | 53.1% | 61.4% |
| `high_confidence` | 14.7% | 41.7% | 51.7% |

<img src="images/benchmark/sensitivity_by_expression_database_human.png" width="100%"/>

Same tertiles, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison:

<img src="images/benchmark/sensitivity_by_expression_database_human_isoform.png" width="100%"/>

#### Isoform-level expression quantification

The table above reports relaxed-BSJ matching, which collapses a locus with several real isoforms into a single point. The tables below instead score each real GT isoform on its own: `full`, `strict`, `relaxed` are the same 3 tiers used for [Isoform detection](#isoform-detection-2) and the isoform PR scatter above (`full` = exact BSJ + 100% exon overlap, `strict` = BSJ within tolerance + 99% overlap, `relaxed` = BSJ within tolerance + 90% overlap). A GT isoform's predicted count is whichever caller isoform structurally matched it under that tier, or 0 if none did, so n here is a genuinely stricter re-match, not the table above's n expanded, and can be smaller, including for the many single-isoform loci in this dataset (isoform identity still has to be structurally confirmed for every one of them, not just the 2 multi-isoform families noted in [Isoform detection](#isoform-detection-2) above). Pearson r (n):

**`full` (exact BSJ, 100% exon overlap):**

| Caller | Database (human) r (n) |
| --- | --- |
| CIRI-long | 0.910 (3,622) |
| isoCirc | 0.553 (5,688) |
| circNICK-lrs | 0.652 (1,470) |
| circFL-seq | 0.776 (10,279) |
| `discovery` | 0.974 (11,369) |
| `balanced` | 0.974 (11,032) |
| `high_confidence` | 0.975 (9,138) |

**`strict` (BSJ within tolerance, 99% exon overlap):**

| Caller | Database (human) r (n) |
| --- | --- |
| CIRI-long | 0.916 (9,081) |
| isoCirc | 0.551 (5,772) |
| circNICK-lrs | 0.660 (1,577) |
| circFL-seq | 0.777 (10,311) |
| `discovery` | 0.974 (12,037) |
| `balanced` | 0.973 (11,697) |
| `high_confidence` | 0.973 (9,252) |

**`relaxed` (BSJ within tolerance, 90% exon overlap):**

| Caller | Database (human) r (n) |
| --- | --- |
| CIRI-long | 0.916 (9,423) |
| isoCirc | 0.551 (5,890) |
| circNICK-lrs | 0.690 (1,870) |
| circFL-seq | 0.777 (10,379) |
| `discovery` | 0.974 (12,352) |
| `balanced` | 0.972 (12,008) |
| `high_confidence` | 0.972 (9,383) |

circNICK-lrs is by far the most affected caller here: its n drops from 5,409 at the locus level above to 1,470-1,870 once every call needs a genuinely matching isoform structure, and its correlation drops with it (r 0.65-0.69 here vs r 0.78 there), unlike the other raw tools and every nanocirc tier, which stay close to the table above at every tier of stringency. CIRI-long's own n also drops hard under `full` specifically (3,622 vs 9,426 at the locus level), less so under `strict`/`relaxed` (9,081-9,423), since exact-structure matching is a much higher bar than boundary tolerance plus 90-99% overlap. Every nanocirc tier still correlates near-perfectly with ground truth (r 0.97-0.98) regardless of stringency, and isoCirc is still the weakest quantifier throughout, the same pattern as every other benchmark in this document.

<img src="images/benchmark/expression_correlation_database_human_full.png" width="100%"/>
<img src="images/benchmark/expression_correlation_database_human_strict.png" width="100%"/>
<img src="images/benchmark/expression_correlation_database_human_relaxed.png" width="100%"/>

<img src="images/benchmark/reproducibility_pairs_database_human_full.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human_strict.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human_relaxed.png" width="100%"/>

### Caller agreement

Ground-truth loci recovered by each of the four tools are pooled across runs and drawn as an UpSet-style plot, where a locus counts as recovered if it was matched in at least one run.

<img src="images/benchmark/upset_plot_database_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/upset_plot_database_human.png" width="100%"/>
<img src="images/benchmark/upset_plot_database_human_exon_based.png" width="100%"/>
<img src="images/benchmark/upset_plot_database_human_full.png" width="100%"/>
<img src="images/benchmark/upset_plot_database_human_strict.png" width="100%"/>
<img src="images/benchmark/upset_plot_database_human_relaxed.png" width="100%"/>

</details>

**What we see:**

- **Database circRNAs are detected far more easily than de novo ones at every tier**: `discovery` reaches F1 0.60 here versus F1 0.47 (human) / 0.56 (mouse) on the synthetic species benchmark, mostly on recall (R 0.48 vs 0.36 / 0.49), consistent with these being already-annotated, pre-selected-for-detectability circRNAs.
- **Quantification is unaffected**: nanocirc tiers still correlate near-perfectly with ground truth (r 0.92-0.95), matching the synthetic benchmarks.
- **Isoform-level numbers here are too small for a trustable result**: only 2 multi-isoform families exist in this database sample, so treat the isoform/family percentages above as descriptive, not something to compare against the other benchmarks' isoform numbers.

**Recommendations**: expect noticeably higher recall (and similar precision) on real, previously-annotated circRNAs than the de novo numbers elsewhere in this document suggest, but keep in mind this dataset is biased toward what tools can already detect (see caveats in the intro above) and don't compare its isoform-level numbers vs other benchmarks given that the randomised database-based sample lacked them.

<details>
<summary><strong>Statistical significance</strong></summary>

### Statistical significance

Do the synthetic benchmark's conclusions transfer to the database-derived dataset? Database circRNAs are easier to detect than synthetic ones by design. They are known, already-catalogued circRNAs, already selected for detectability, and used as input or a built-in dependency by tools during circRNA calling. The in silico dataset instead tries to evaluate performance on novel, unannotated molecules: it asks "if the read evidence for this molecule is there, but no annotation hints exist, can the tools still recover it?" The database benchmark asks a different question instead: "of circRNAs that are already known, annotated, and previously detected by different tools, how many can a tool recover when given read evidence and, where possible, the same input database?"

So the gap between the two methods is expected, but it's hard to say how much of it is due to the simulation logic, the database being biased by pre-selection for detectability, the presence of circRNA structures that were already reported by the same tools (which may still be real TPs or FPs, since we don't have wet-lab validation data for each database entry and depositing one doesn't require it), the presence of a database input file that directly hints the tools, actual biology, or other factors. What we can measure is whether tool ranking and per-type detectability can be transferred between the two. 

We matched precision/recall/F1 between the in silico and database benchmarks at strict isoform-level matching (BSJ ±10bp and ≥99% exon overlap), separately for every caller, circRNA type, and run: up to 120 paired points (8 callers × 5 types × 3 runs; precision has 111, since a caller that made zero calls of a type in a run has undefined precision, not 0, so that cell is dropped). Computed by [`scripts/statistical_analysis.py`](https://gitlab.com/bingo-igdr/nano-circ), reusable for any two datasets scored with our benchmark (not just in silico vs database, also wet-lab protocol vs protocol and species vs species, both used elsewhere in this document).

Same statistical methods as the [wet-lab-protocol comparison](#statistical-significance) (linear regression, Kendall's tau-b, Mantel-Haenszel + Breslow-Day, Wilcoxon signed-rank, Cohen's d), applied here to synthetic vs database instead of protocol vs protocol.

#### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | Synthetic | Database | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 99 | 0.705 | 0.830 | +0.125 | 0.696 | 0.339 | 0.607 | 4×10⁻¹⁶ | 0.726 |
| Recall | 105 | 0.241 | 0.324 | +0.084 | 0.745 | 0.145 | 0.644 | 2×10⁻⁹ | 0.651 |
| F1 | 105 | 0.305 | 0.435 | +0.130 | 0.810 | 0.188 | 0.584 | 8×10⁻¹⁴ | 0.872 |

Recall and F1 on database circRNAs are well predicted by their in silico counterpart values (R² 0.58–0.64), a shrinking-toward-the-middle relationship (slope 0.75–0.81, below 1): callers that were already strong on synthetic circRNAs gain proportionally less than weaker ones. Precision is similarly predictable (R² = 0.61, slope = 0.70). All three metrics are significantly higher on database circRNAs (Wilcoxon p < 10⁻⁸, Cohen's d 0.65–0.87).

#### Per type, by metric (n=21 paired points per row pooling all 7 callers × 3 runs, except antisense precision, n=21, and intergenic precision, n=18: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | Synthetic | Database | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 21 | 0.629 | 0.813 | +0.184 | 0.899 | 0.248 | 0.843 | 1×10⁻⁶ | 2.169 |
| eciRNA | Recall | 21 | 0.515 | 0.558 | +0.043 | 1.097 | −0.007 | 0.933 | 9×10⁻⁴ | 1.006 |
| eciRNA | F1 | 21 | 0.553 | 0.656 | +0.104 | 1.104 | 0.046 | 0.952 | 1×10⁻⁶ | 2.566 |
| EIciRNA | Precision | 21 | 0.650 | 0.762 | +0.112 | 1.118 | 0.035 | 0.872 | 7×10⁻⁴ | 0.984 |
| EIciRNA | Recall | 21 | 0.113 | 0.149 | +0.036 | 0.975 | 0.039 | 0.698 | 4×10⁻⁴ | 1.042 |
| EIciRNA | F1 | 21 | 0.189 | 0.248 | +0.059 | 1.035 | 0.052 | 0.748 | 4×10⁻⁴ | 1.142 |
| ciRNA | Precision | 18 | 0.882 | 0.947 | +0.065 | 0.361 | 0.629 | 0.857 | 1×10⁻⁴ | 1.161 |
| ciRNA | Recall | 21 | 0.118 | 0.328 | +0.210 | 1.795 | 0.117 | 0.648 | 1×10⁻⁶ | 1.893 |
| ciRNA | F1 | 21 | 0.198 | 0.461 | +0.264 | 1.628 | 0.139 | 0.740 | 1×10⁻⁶ | 2.034 |
| antisense | Precision | 21 | 0.648 | 0.836 | +0.188 | 0.490 | 0.519 | 0.410 | 1×10⁻⁶ | 0.580 |
| antisense | Recall | 21 | 0.157 | 0.318 | +0.161 | 0.583 | 0.226 | 0.256 | 1×10⁻⁶ | 1.156 |
| antisense | F1 | 21 | 0.245 | 0.455 | +0.210 | 0.539 | 0.323 | 0.254 | 1×10⁻⁶ | 1.068 |
| intergenic | Precision | 18 | 0.746 | 0.803 | +0.057 | 1.012 | 0.048 | 0.933 | 3×10⁻⁴ | 1.226 |
| intergenic | Recall | 21 | 0.300 | 0.269 | −0.030 | 0.656 | 0.073 | 0.904 | 0.215 | −0.292 |
| intergenic | F1 | 21 | 0.341 | 0.355 | +0.014 | 0.728 | 0.106 | 0.826 | 0.744 | 0.135 |

Every type's recall and F1 transfer well between the two benchmarks (R² 0.25–0.95): a caller that's relatively better or worse at a type on synthetic circRNAs stays that way on database circRNAs. intergenic is the one type where recall/F1 show no significant database boost at all (p = 0.22 and 0.74), the only type in this table where de novo and database performance are statistically indistinguishable on those two metrics. Every other type/metric combination is significantly higher on database circRNAs (p ≤ 7×10⁻⁴, Cohen's d 0.58–2.57), including every type's precision (p ≤ 1×10⁻⁴).

#### Does caller ranking transfer?

We averaged each caller's synthetic and database values across type and run first, ranked the 7 callers by that average on each side, then compared the two rankings with Kendall's tau-b:

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 7 | 0.619 (p = 0.069) | 0.750 (p = 0.052) |
| Recall | 7 | 0.714 (p = 0.030) | 0.857 (p = 0.014) |
| F1 | 7 | 0.714 (p = 0.030) | 0.857 (p = 0.014) |

Recall/F1 ranking transfers significantly (tau-b = 0.71); precision ranking is weaker and not quite significant (tau-b = 0.62, p = 0.069). By F1:

| Caller | F1 synthetic | F1 database | Rank synthetic | Rank database |
| --- | --- | --- | --- | --- |
| balanced | 0.445 | 0.581 | 1 | 1 |
| discovery | 0.428 | 0.563 | 2 | 2 |
| CIRI-long | 0.375 | 0.470 | 3 | 5 |
| high_confidence | 0.329 | 0.491 | 4 | 4 |
| circFL-seq | 0.302 | 0.539 | 5 | 3 |
| isoCirc | 0.206 | 0.338 | 6 | 6 |
| circNICK-lrs | 0.050 | 0.063 | 7 | 7 |

Five of the seven callers had the exact same rank on both benchmarks. circFL-seq and CIRI-long swapped places: circFL-seq ranked 5th on in silico circRNAs but 3rd on database ones, CIRI-long the reverse. This suggests that circFL-seq is stronger on known annotated circRNAs and worse on de novo.

<img src="images/benchmark/slopegraph_database_human.png" width="100%"/>

#### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=9)

| Tool | Metric | n | Synthetic | Database | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.695 | 0.865 | +0.170 | 0.214 | 0.717 | 0.108 | 6×10⁻⁵ | 1.514 |
| CIRI-long | Recall | 15 | 0.293 | 0.341 | +0.047 | 0.785 | 0.111 | 0.737 | 0.421 | 0.493 |
| CIRI-long | F1 | 15 | 0.375 | 0.470 | +0.095 | 1.169 | 0.031 | 0.681 | 0.026 | 0.934 |
| isoCirc | Precision | 15 | 0.886 | 0.942 | +0.056 | 0.420 | 0.570 | 0.900 | 0.003 | 0.920 |
| isoCirc | Recall | 15 | 0.135 | 0.216 | +0.081 | 0.717 | 0.119 | 0.888 | 6×10⁻⁵ | 1.445 |
| isoCirc | F1 | 15 | 0.206 | 0.338 | +0.133 | 0.627 | 0.209 | 0.840 | 6×10⁻⁵ | 1.450 |
| circNICK-lrs | Precision | 9 | 0.091 | 0.146 | +0.055 | 1.308 | 0.027 | 0.554 | 0.055 | 0.605 |
| circNICK-lrs | Recall | 15 | 0.057 | 0.056 | −0.001 | 1.011 | −0.002 | 0.998 | 0.875 | −0.197 |
| circNICK-lrs | F1 | 15 | 0.050 | 0.063 | +0.013 | 1.382 | −0.006 | 0.982 | 0.239 | 0.370 |
| circFL-seq | Precision | 15 | 0.679 | 0.924 | +0.245 | −0.053 | 0.960 | 0.076 | 6×10⁻⁵ | 0.660 |
| circFL-seq | Recall | 15 | 0.212 | 0.393 | +0.180 | 0.452 | 0.297 | 0.357 | 6×10⁻⁵ | 1.186 |
| circFL-seq | F1 | 15 | 0.302 | 0.539 | +0.237 | 0.261 | 0.460 | 0.154 | 6×10⁻⁵ | 1.139 |
| high_confidence | Precision | 15 | 0.900 | 0.960 | +0.059 | 0.369 | 0.627 | 0.556 | 6×10⁻⁵ | 1.115 |
| high_confidence | Recall | 15 | 0.231 | 0.347 | +0.116 | 0.583 | 0.212 | 0.438 | 0.015 | 0.827 |
| high_confidence | F1 | 15 | 0.329 | 0.491 | +0.161 | 0.411 | 0.355 | 0.242 | 0.008 | 0.884 |
| balanced | Precision | 15 | 0.741 | 0.886 | +0.145 | 0.540 | 0.486 | 0.603 | 6×10⁻⁵ | 2.077 |
| balanced | Recall | 15 | 0.374 | 0.454 | +0.080 | 0.524 | 0.258 | 0.537 | 0.135 | 0.517 |
| balanced | F1 | 15 | 0.445 | 0.581 | +0.136 | 0.542 | 0.340 | 0.440 | 9×10⁻⁴ | 0.992 |
| discovery | Precision | 15 | 0.696 | 0.811 | +0.114 | 0.957 | 0.144 | 0.714 | 9×10⁻⁴ | 1.216 |
| discovery | Recall | 15 | 0.381 | 0.465 | +0.084 | 0.533 | 0.262 | 0.561 | 0.135 | 0.542 |
| discovery | F1 | 15 | 0.428 | 0.563 | +0.135 | 0.540 | 0.332 | 0.343 | 9×10⁻⁴ | 0.942 |

`high_confidence` has full n=15 coverage on precision, consistent with the wet-lab comparisons above. isoCirc and `high_confidence` are consistent across all three metrics (R² ≥ 0.55). `balanced` and `discovery` are strong on precision (R² ≥ 0.60) but less consistent on recall and F1 (R² 0.34–0.56, not always significant). CIRI-long and circNICK-lrs are the weakest fits: several of their rows are not significant at all (CIRI-long's recall p = 0.42, circNICK-lrs's recall/F1 p ≥ 0.24). By F1, circNICK-lrs and CIRI-long gain the least on database circRNAs (+0.01 to +0.10), while `balanced`, `discovery`, and circFL-seq gain the most (+0.14 to +0.24).

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 3.72 | 3.47–3.99 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| CIRI-long | Recall | 5 | 1.29 | 1.24–1.34 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 2.13 | 1.86–2.43 | <10⁻³⁰⁰ | 0.022 |
| isoCirc | Recall | 5 | 1.83 | 1.74–1.92 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circNICK-lrs | Precision | 5 | 1.84 | 1.70–1.99 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circNICK-lrs | Recall | 5 | 0.98 | 0.90–1.06 | 0.586 | 3×10⁻⁵ |
| circFL-seq | Precision | 5 | 2.67 | 2.43–2.92 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circFL-seq | Recall | 5 | 2.59 | 2.49–2.70 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| high_confidence | Precision | 5 | 3.05 | 2.71–3.42 | <10⁻³⁰⁰ | 3×10⁻¹⁶ |
| high_confidence | Recall | 5 | 1.92 | 1.85–2.00 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced | Precision | 5 | 2.84 | 2.67–3.03 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced | Recall | 5 | 1.48 | 1.43–1.54 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 1.78 | 1.69–1.87 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Recall | 5 | 1.52 | 1.47–1.58 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

Every caller has significantly higher odds of a correct call (precision) and of recovering a ground-truth circRNA (recall) on database circRNAs than on synthetic ones (OR 1.3–3.7), except circNICK-lrs's recall, which shows no real difference (OR ≈ 0.98, matching its flat recall line in the slopegraph). But that boost isn't the same size across circRNA types: of the fourteen caller/outcome combinations tested, none are uniform across types (Breslow-Day p ≤ 0.022 for every row); every combination's boost depends on which circRNA type it is.

**Conclusion.** Correlation is reasonable: aggregate R² is 0.58–0.64, and most types and callers individually fit a significant line, though CIRI-long's and circNICK-lrs's per-tool fits are weaker. Caller ranking transfers on recall/F1 (Kendall's tau-b = 0.71) but not quite on precision (tau-b = 0.62, p = 0.069); five of the seven callers land on the same F1 rank on both benchmarks, with only circFL-seq and CIRI-long swapping places with each other. The size of the database boost isn't uniform: it depends on circRNA type for every one of the fourteen caller/outcome combinations tested. So the overall pattern and ranking are preserved reasonably well between de novo and database calling, but the exact performance improvement on annotated circRNAs is not uniform, and depends on both the circRNA type and the caller.

</details>

<details>
<summary><strong>What are these false positives?</strong></summary>

### What are these false positives?

We traced every FP-supporting read back to the simulator's own ground truth (single-tool discovery exon-based FPs, all 3 human runs, mean % ± SD across runs):

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 29.0 ± 1.0 | 64.8 ± 2.1 | 88.8 ± 2.5 | 0.6 ± 0.1 |
| Real circ read, correct strand, pure boundary noise | 20.8 ± 2.6 | 16.5 ± 2.5 | 0.0 ± 0.0 | 34.4 ± 0.6 |
| Real circ read, but circRNA reported at a different/distant locus | 50.3 ± 1.9 | 18.7 ± 0.4 | 11.2 ± 2.5 | 65.0 ± 0.6 |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.1 ± 0.1 |
| FP locus overlaps a genomic repeat | 57.7 ± 1.7 | 56.4 ± 1.1 | 47.4 ± 1.4 | 87.9 ± 0.3 |
| n FP loci (pooled, 3 runs) | 635 | 653 | 118 | 4,930 |
| n supporting reads (pooled, 3 runs) | 12,614 | 5,608 | 1,052 | 47,269 |

This is a genuinely different FP composition from every synthetic-data comparison above. On synthetic data, wrong-strand antisense calls dominate every tool's FPs (77-94%, see the wet-lab-protocol and species sections' equivalent tables). On database circRNAs, only isoCirc still shows that pattern (88.8%); CIRI-long and circNICK-lrs are instead dominated by different-locus calls (50.3% and 65.0%), and circFL-seq is roughly split between wrong-strand (64.8%) and boundary/locus noise. Repeat overlap is still high across all four tools (47-88%), as it is throughout this document. A plausible explanation is that unlike the simulator, the database annotation carries circRNAs with more diverse, sometimes ambiguous, genomic contexts (isoforms of overlapping/nearby genes, non-simulated antisense pairs), so a caller can land near the right region without hitting the exact same locus or transcribing from the correct strand's known circRNA.

**Practical implication**: database `discovery` false positives are less concentrated in a single, explainable error mode than synthetic ones. Re-verify low-confidence calls the same way recommended elsewhere in this document (e.g. Sanger sequencing) before trusting exact structure, strand, or locus.

</details>

<details>
<summary><strong>What discovery still misses?</strong></summary>

### What discovery still misses?

We also characterized every ground-truth circRNA against whether any tool recovered it at all (mean % ± SD across 3 human runs; missed + recovered = 100% of the GT catalog, ~8,672 circRNAs/run):

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 52.1 ± 0.1 | 47.9 ± 0.1 |
| Mean mature length (bp) | 13,292 ± 66 | 319 ± 0 |
| Mean rolling-circle copy number | 7.73 ± 0.05 | 8.54 ± 0.04 |
| Median simulated read count | 31 ± 0 | 82 ± 1 |
| Canonical GT/AG splice site at BSJ | 94.8 ± 0.1 | 96.3 ± 0.0 |
| No gene overlap at all (intergenic) | 16.4 ± 0.5 | 15.9 ± 0.4 |
| Falls in a gene with ≥5 isoforms | 72.2 ± 0.6 | 71.9 ± 0.4 |
| Overlaps a genomic repeat | 78.5 ± 0.0 | 78.8 ± 0.1 |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 25.9 ± 0.2 | 74.1 ± 0.2 |
| `ciRNA` (purely intronic) | 49.5 ± 0.4 | 50.5 ± 0.4 |
| `antisense` | 51.1 ± 0.1 | 48.9 ± 0.1 |
| `intergenic` | 64.7 ± 0.1 | 35.3 ± 0.1 |
| `EIciRNA` (exon+intron) | 72.8 ± 0.0 | 27.2 ± 0.0 |

Missed circRNAs are far longer (~42x, versus 4-9x on the synthetic benchmarks) and lower-expressed, but the copy-number, splice-site, intergenic, and repeat-overlap gaps between missed and recovered are all small here, much smaller than the synthetic benchmarks' equivalents. `eciRNA` is again the easiest biotype; `EIciRNA` the hardest, matching the type-classification table above and every other benchmark in this document.

</details>

---

## Cross-benchmark comparison: which effect is biggest?

Every benchmark above shows two datasets with difference in one variable (species, wet-lab protocol, or annotated-vs-de-novo). This section plots them in a single PCA over each dataset's pooled precision/recall/F1, computed by [`publication/pca_analysis.py`](https://gitlab.com/bingo-igdr/nano-circ): one row per dataset, one column per (caller × structural strategy × threshold × metric) combination shared by all 7, standardized (z-scored) before decomposition.

<img src="images/benchmark/pca_benchmarks.png" width="100%"/>

PC1 (72% of variance) is dominated by F1 and recall, i.e. overall detection difficulty (how hard to make a call); PC2 (15%) is dominated by precision (how precise are made calls).

To isolate each effect on its own, `pca_analysis.py` also reports the direct distance (in the same standardized feature space the PCA was computed in, not just the PC1/PC2 projection) between pairs of datasets differing in only that one factor, saved to `pca_pairwise_distances.csv`:

| Effect | Datasets compared | Distance |
| --- | --- | --- |
| Annotated (database) vs de novo (synthetic) | database_human vs species_human | 19.13 |
| isoCirc vs circFL-seq protocol | wetlab_isocirc vs wetlab_circfl | 15.70 |
| Human vs mouse | species_human vs species_mouse | 15.48 |
| isoCirc vs CIRI-long protocol | wetlab_isocirc vs wetlab_ciri_long | 12.75 |
| CIRI-long vs circFL-seq protocol | wetlab_ciri_long vs wetlab_circfl | 6.75 |

**Conclusion.** Among everything tested, the biggest variation comes from annotated (database) vs de novo (synthetic) circRNA recovery, with annotated circRNAs being easier to identify, closely followed by the isoCirc-protocol wet-lab effect. The species effect (human vs mouse) is smaller but in the same range as the isoCirc-protocol effect. CIRI-long vs circFL-seq protocol has the smallest effect of anything tested, under half the size of every other comparison. 

The annotated vs de novo distance (≈19) is not disproportionately large next to other, it's in the same order of magnitude. That supports the earlier hypothesis that de novo synthetic circRNAs are a valid comparator: they're offset from database circRNAs in absolute difficulty (expected as per reasons already documented above), by a margin similar in scale to other differences (species, protocol) while caller ranking and per-type structure transfer well between the two (see [Statistical significance](#statistical-significance-2) above).

Putting all benchmarks together: annotated circRNAs are both more easily recoverable and called with higher precision than de novo ones. Among wet-lab protocols, circFL-seq allows the best recovery, closely followed by CIRI-long protocol; isoCirc protocol is the hardest for the tools to recover reads from, despite offering slightly higher precision. Species differences are also to be expected: mouse-derived circRNAs were easier to recover than human-derived ones, with an effect similar in size to a wet-lab-protocol effect (species distance 15.48, close to isoCirc-vs-circFL-seq's 15.70). Tool ranking generally transfers well across every one of these comparisons (Kendall's tau-b 0.71–1.00), and when a discrepancy exists it's usually just two callers swapping adjacent places, for example circFL-seq and CIRI-long swapping on annotated circRNAs, or circFL-seq and isoCirc swapping between human and mouse. Across all of this, nanocirc's own tiers are the most stable callers: `balanced` and `discovery` held the top two F1 ranks in every wet-lab-protocol, species, and database comparison in this benchmark, and every nanocirc tier correlated with ground-truth expression more strongly than any single raw tool on every dataset tested, while the raw tools' relative strengths shifted more with protocol and species.

---

For output file formats see [docs/output.md](output.md). For parameters see [docs/usage.md](usage.md). For merge/scoring algorithms see [docs/methods.md](methods.md).
