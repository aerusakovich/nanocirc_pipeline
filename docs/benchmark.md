# nanocirc: Benchmark Results and Recommendations

This benchmark is adapted from our own earlier benchmark study, Rusakovich et al., "Benchmarking circRNA detection tools from long-read sequencing using a data-driven and flexible simulation framework" (Peer Community Journal, 2026; see [CITATIONS.md](../CITATIONS.md) for the full reference), modified and extended based on our newer findings. 

This page summarizes the benchmark behind nanocirc's four confidence tiers and gives practical guidance on which one to use. The new simulator and benchmark scripts is in this repo, [nano-circ](https://gitlab.com/bingo-igdr/nano-circ), as v3. Results below are from simulated human and mouse datasets, 3 runs each.

Several benchmarks are shown below:

- **[Same circRNAs, three wet-lab protocols](#same-circrnas-three-wet-lab-protocols)**: NanoSim's read-error model is characterized from real CIRI-long/isoCIRC and circFL-seq protocols ONT control reads; circRNA type composition and length/exon/splice distributions are the same across all 3 wet-lab methods, circRNA generation based on circRNAs from a human circbase+circatlas catalog. 
- **[Species benchmark: human vs mouse](#species-benchmark-human-vs-mouse)**: CIRI-long protocol, each species' own circRNA catalog and genome. Any difference here comes from species biology (genome, catalog composition), not the wet-lab protocol, which is fixed to CIRI-long for both.
- **[Database benchmark: annotated circRNA recovery](#database-benchmark-annotated-circrna-recovery)**: annotated vs de novo circRNA recovery. Same simulation and pipeline as above, but circRNAs are sampled from a real human circBase/circAtlas catalog instead of generated de novo, so tools are evaluated on already-known, already-catalogued circRNAs instead of novel unannotated ones.

## Recommendations

- **Choose circFL-seq wet-lab protocol if you can, followed by CIRI-long, for maximum tool compatibility.** Both work with all 4 tools; circFL-seq protocol gave the highest recall, highest F1, and the highest BSJ-confirmation rate of the 3 protocols we tested (64.1% vs CIRI-long's 54.8%), and its read structure also allowed to better resolve isoform structure. isoCirc-protocol reads are the hardest for the tools to work with, and CIRI-long specifically is almost incompatible with them. See [Structural precision/recall](#structural-precisionrecall) and [Isoform detection](#isoform-detection) under [Same circRNAs, three wet-lab protocols](#same-circrnas-three-wet-lab-protocols).
- **Use nanocirc itself rather than any single tool, and use its own quantification (`--run_quantify`) rather than a raw tool's built-in read counting.** All 4 tools are containerized and patched against real bugs found running them at scale, several silent rather than crashing (wrong fastq/genome argument wiring, a BAM filename mismatch that silently dropped intron-coverage output, a BioPython-version crash that only crashes on certain read alignments). A tool being good at *finding* the right circRNA doesn't mean its own built-in quantifier counts it well and native quantifiers can't quantify circRNAs tool didn't report. nanocirc's remap-based quantification is not dependent on which tool/tier found a locus. It applies the same counting method regardless of source, is stable across wet-lab approaches and has been benchmarked to agree well across tools/tiers and against ground truth (see [Quantification advice](#quantification-advice)). It also builds one shared, cross-sample circRNA catalog (`--run_crossrun_merge`) for a multi-sample count matrix that can be used with DESeq2 directly.
- **For isoform-level differential expression, use the `balanced_recall` tier.** It is close to `discovery`'s isoform-family resolution while keeping meaningfully higher precision, so you get `discovery`'s isoform sensitivity with fewer false positives in the count matrix. `discovery` is an acceptable choice for locus-level differential expression, where getting the exact internal exon structure right isn't the point. See [Isoform detection](#isoform-detection).
- **Exception: on isoCirc-protocol data, use `discovery` for isoform-level work instead of `balanced_recall`.** circNICK-lrs's own structure calls aren't trustworthy enough to count as consensus evidence, and CIRI-long is close to incompatible with isoCirc-protocol reads (see [Which tools work with my wet-lab protocol and species?](#which-tools-work-with-my-wet-lab-protocol-and-species)), so the only tools left to agree on a recovered isoform are circFL-seq and isoCirc, which rarely call the exact same structure independently. In our isoCirc-protocol benchmarks, `balanced_recall`'s isoform-family resolution dropped to 0%, while `discovery` still reached ~12%. See [Isoform detection](#isoform-detection).
- **To verify a specific circRNA's exon boundaries before follow-up work (primers, probes, isoform-specific validation), check it against the `balanced_precision`/`high_confidence` sets.** Both only contain calls with strong cross-tool agreement on both BSJ and structure, so an entry present there has cleared more quality thresholds. See [Which tier should I use?](#which-tier-should-i-use) and [Isoform detection](#isoform-detection).

## Which tools work with my wet-lab protocol and species?

nanocirc runs 4 tools (CIRI-long, circFL-seq, isoCirc, circNICK-lrs) by default. Not every tool works with every wet-lab protocol's read structure, and circNICK-lrs specifically is also species-dependent, not just protocol-dependent (see table). Set `--wet_lab` to your protocol (`ciri_long`, `circfl`, `isocirc`, or `circnick`) to use default settings for this type of data. `--wet_lab circnick` switches to circNICK-lrs alone: this is the one protocol where the other tools would certainly fail, since the other 3 tools expect rolling-circle reads with multiple copies of the circRNA per read, which is not the case for circNICK-lrs wet-lab data. `--wet_lab isocirc` turns CIRI-long off by default (see below). For every other protocol, `--wet_lab` does not turn any tool off by itself: so we suggest you to try with the default settings, but you can turn on/off tools with `--run_isocirc`, `--run_circfl`, `--run_circnick`, `--run_cirilong` set to `false`.

| Your wet-lab protocol | Tools suggestion |
| --- | --- |
| circNICK-lrs | circNICK-lrs only, human/mouse only**** |
| CIRI-long | All 4 tools, circNICK-lrs only if human/mouse**** |
| circFL-seq | All 4 tools, circNICK-lrs only if human/mouse**** |
| isoCirc* | 3 tools, CIRI-long off by default**, circNICK-lrs only if human/mouse**** |
| PacBio | protocol-dependent, turn CIRI-long off manually if isoCirc based***, circNICK-lrs only if human/mouse**** |


\* circFL-seq also needs more memory than usual on isoCirc data: a first attempt at 250GB ran out of memory, a second attempt completed at 500GB, after about 5.5 days. Give circFL-seq extra memory if your samples come from isoCirc-protocol reads.

\*\* On real isoCirc-protocol CIRI-long did not finish in reasonable time, getting stuck on its collapse step for days. We tried it 3 times (5 days, then 10 days, then over 24 days), it still has not completed. It worked on simulated data and 10 pct of the real isoCirc dataset, but much slower than on the datasets from other wet-lab protocols, so this issue may be dataset size dependent. To try CIRI-long on isoCirc data use `--run_cirilong true`. By our benchmark, leaving it off on isoCirc data costs about 0.01 F1 at the `discovery`, `balanced_recall`, and `balanced_precision` tiers, and actually improves `high_confidence` (F1 almost doubles) while improving compute time and reducing the risk of pipeline freeze on this tool. 

\*\*\* On real isoCirc-based PacBio human data, CIRI-long did not finish in reasonable time, getting stuck on its collapse step for days. We tried it 3 times (5 days, then 10 days, then over 24 days), it still has not completed. It worked on simulated data and 10 pct of the real isoCirc-based PacBio dataset, so this issue may be dataset size dependent. PacBio has no dedicated `--wet_lab` preset, as it is a sequencing platform, so select the preset based on the wet-lab approach and turn `--run_cirilong` off manually if you hit the same issue.

\*\*\*\* circNICK-lrs is only benchmarked and supported for human and mouse. Turn it off with `--run_circnick false` for any other species, we have not validated it elsewhere and it is not expected to work.

## Which tier should I use?

Each tier trades recall against precision differently. There is no single best tier, only the right one for your question.

| Tier | What it optimizes for | Trade-off |
| ---- | ---------------------- | --------- |
| **`discovery`** | Recall. Keeps almost every circRNA any tool found. | Highest false-positive rate: computational artefacts and incorrectly predicted structures will pass alongside real calls. |
| **`high_confidence`** | Precision. Only calls with strong cross-tool agreement on both BSJ and structure or isoCirc calls that met defined read threshold survive. | Sacrifices recall to get the structure right: false positives are minimal, but the reduced set will miss real biology so it's not fit for exploratory studies. |
| **`balanced_recall`** | A middle ground leaning recall. Cuts most of `discovery`'s false-positive rate while keeping reasonable recall, including on isoform structure (see [Isoform detection](#isoform-detection)). | Good default for exploratory analyses and for isoform-level questions, where `discovery` alone would carry more false positives for similar isoform resolution. |
| **`balanced_precision`** | A middle ground leaning precision. Keeps recall relatively high while pushing precision further than `balanced_recall`. | Best for precision-leaning analyses that need to avoid false positives but still want decent recall, not just the aggressively small `high_confidence` set. |

> **Recommendation:** use one of the two **balanced** tiers for most applications. Reserve `discovery` for exploratory, locus-level studies, and `high_confidence` for when you specifically need a small, near-certain set.

## Quantification advice

Quantification is its own remap step, independent of which tool originally called a locus: the reproducibility and expression correlation plots below show read counts agree well across tools and tiers, and against ground truth, so switching tier for a re-analysis doesn't mean changing to a different counting method.

Please note that for differential expression (read counts, DESeq2), `discovery`'s false positives are not easy-to-spot near-zero-read artefacts that low-count filtering could remove (median ~16 reads vs ~16-23 for true positives in our benchmark), so they will be present in the count matrix as real-looking features. `high_confidence`, on the other end, drops enough real calls to hurt power on lower-expressed circRNAs. This is why we recommend the Balanced modes for expression-based analyses.

For a multi-sample design, turn on `--run_crossrun_merge` alongside `--run_quantify` so every sample in a `group` is quantified against one shared catalog. Otherwise each sample gets its own catalog and rows don't line up across samples in the [DESeq2 matrix](output.md#deseq2-count-matrix).

## Same circRNAs, three wet-lab protocols

We used one circRNA catalog (same circbase/circatlas composition, human genome, same abundances) and only changed the wet-lab protocol NanoSim is characterized from. Any difference in the results below comes from the wet-lab protocol characteristics.

<details>
<summary><strong>Structural precision/recall</strong></summary>

### Structural precision/recall

<img src="images/benchmark/pr_scatter_equal_ciri_long.png" width="100%"/>
<img src="images/benchmark/pr_scatter_equal_circfl.png" width="100%"/>
<img src="images/benchmark/pr_scatter_equal_isocirc.png" width="100%"/>

Same circRNA catalog, same tools, same tiers, only the protocol the reads were simulated from changes between these 3 plots.

Mean precision (P) / recall (R) / F1 across 3 runs, exon-based matching, per nanocirc tier, plus each individual tool's own raw calls for comparison:

| Tier / tool | CIRI-long protocol | circFL-seq protocol | isoCirc protocol (CIRI-long on) |
| --- | --- | --- | --- |
| isoCirc alone | P 0.89 / R 0.14 / F1 0.25 | P 0.90 / R 0.19 / F1 0.32 | P 0.87 / R 0.08 / F1 0.15 |
| circFL-seq alone | P 0.79 / R 0.21 / F1 0.33 | P 0.77 / R 0.25 / F1 0.37 | P 0.74 / R 0.23 / F1 0.36 |
| CIRI-long alone | P 0.74 / R 0.34 / F1 0.46 | P 0.77 / R 0.41 / F1 0.53 | P 0.77 / R 0.03 / F1 0.05 |
| circNICK-lrs alone | P 0.24 / R 0.09 / F1 0.13 | P 0.24 / R 0.10 / F1 0.14 | P 0.22 / R 0.10 / F1 0.14 |
| `discovery` | P 0.68 / R 0.43 / F1 0.53 | P 0.69 / R 0.52 / F1 0.59 | P 0.66 / R 0.26 / F1 0.37 |
| `balanced_recall` | P 0.74 / R 0.42 / F1 0.53 | P 0.77 / R 0.51 / F1 0.61 | P 0.81 / R 0.25 / F1 0.38 |
| `balanced_precision` | P 0.85 / R 0.22 / F1 0.35 | P 0.85 / R 0.28 / F1 0.42 | P 0.85 / R 0.09 / F1 0.17 |
| `high_confidence` | P 0.86 / R 0.07 / F1 0.14 | P 0.88 / R 0.10 / F1 0.17 | P 0.87 / R 0.01 / F1 0.02 |

</details>

<details>
<summary><strong>Ground truth coverage</strong></summary>

### Ground truth coverage

<img src="images/benchmark/equal_gt_coverage_stacked_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_gt_coverage_stacked_circfl.png" width="100%"/>
<img src="images/benchmark/equal_gt_coverage_stacked_isocirc.png" width="100%"/>

circFL-seq protocol has led to the highest ground truth recovery at every tier (e.g. `discovery` R 0.52 vs CIRI-long's 0.43 and isoCirc's 0.26, from the table above); isoCirc protocol has the lowest at every tier, meaning it's reads are hardest to recover the calls from.

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

Same TP/FP-vs-GT breakdown as the GT coverage plot above, split by circRNA biotype. eciRNA is the easiest type to recover on every protocol. circFL-seq protocol leads to the best recovery of every type; isoCirc protocol to the worst recovery, most severely affecting antisense circRNA (8.4% vs CIRI-long's 32.7% and circFL-seq's 37.8%).

</details>

<details>
<summary><strong>Isoform detection</strong></summary>

### Isoform detection

The structural plots above show if a specific part of circRNA structure was found at its locus. But we are interested to see if the callers have found the exactly same molecule, especially on isoform level, which is the main point of using long-read sequencing. We evaluated that with three tiers, each pairing a BSJ boundary check with an exon overlap check: `full` (exact BSJ, 100% overlap, byte identical structure), `strict` (BSJ within tolerance window, 99% exon overlap, allows little boundary noise but can be easily reconstructed by follow up analyses), `relaxed` (BSJ within tolerance, 90% overlap, allowing a part of mismatched sequence).

<img src="images/benchmark/pr_scatter_isoform_ciri_long.png" width="100%"/>
<img src="images/benchmark/pr_scatter_isoform_circfl.png" width="100%"/>
<img src="images/benchmark/pr_scatter_isoform_isocirc.png" width="100%"/>

Long-read protocols are used mainly to resolve internal exon structure, not just find BSJ position or general locus. Among ground-truth loci where the same BSJ has more than one real isoform (isoform family), we checked how often each caller recovers the exact isoform with different levels of stringency.

<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_99_isocirc.png" width="100%"/>

<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_detection_reciprocal_overlap_90_isocirc.png" width="100%"/>

Green+grey sum to 100% of the isoform-family ground truth; red shows the calls that found the right BSJ but didn't match required exon structure threshold. 99% reciprocal exon overlap allows only trivial boundary noise; 90% allows a small sequence mismatch while still requiring generally the same structure, not a different one. 

% of ground-truth isoforms called at 99% reciprocal overlap:

| Tier | CIRI-long | circFL-seq | isoCirc |
| --- | --- | --- | --- |
| `discovery` | 61.8% | 67.6% | 41.2% |
| `balanced_recall` | 55.9% | 61.8% | 29.4% |
| `balanced_precision` | 38.2% | 44.1% | 23.5% |
| `high_confidence` | 29.4% | 23.5% | 0.0% |

circFL-seq protocol allows resolution of the most isoforms at every tier except `high_confidence`, where CIRI-long is actually ahead (29.4% vs 23.5%). isoCirc protocol is last everywhere, and `high_confidence` never gets a single exact isoform right on it. Relaxing to 90% overlap barely moves these numbers (0-6 percentage points difference).

<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_99_isocirc.png" width="100%"/>

<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_circfl.png" width="100%"/>
<img src="images/benchmark/equal_isoform_family_completeness_reciprocal_overlap_90_isocirc.png" width="100%"/>

Per family plot instead of per isoform: did the caller get every isoform in the family right, some of them, or none. % of isoform families fully resolved (every isoform in the family called exactly right), at 99% reciprocal overlap:

| Tier | CIRI-long | circFL-seq | isoCirc |
| --- | --- | --- | --- |
| `discovery` | 35.3% | 41.2% | 11.8% |
| `balanced_recall` | 23.5% | 29.4% | 0.0% |
| `balanced_precision` | 0.0% | 5.9% | 0.0% |
| `high_confidence` | 0.0% | 0.0% | 0.0% |

Getting one isoform right while missing or misclassifying others is the most common outcome; fully resolving an isoform family is rare for every tier and protocol. circFL-seq protocol gives the best odds of a full family at every tier; isoCirc protocol remains the hardest for all tools.

</details>

<details>
<summary><strong>Expression quantification</strong></summary>

### Expression quantification

<img src="images/benchmark/equal_expression_correlation_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_circfl.png" width="100%"/>
<img src="images/benchmark/equal_expression_correlation_isocirc.png" width="100%"/>

Tool may call a correct circRNA, but misquantify it. It is important to measure, as for many research tasks that rely on differential expression precise circRNA quantification on locus and isoform is required. We measured it with Pearson r, 
Spearman rho, R², read-count correlation against ground truth and cross-method agreement, per tool/tier. 

Pearson r by protocol (relaxed-BSJ matching):

| Caller | CIRI-long | circFL-seq | isoCirc |
| --- | --- | --- | --- |
| CIRI-long | 0.911 | 0.935 | not detected |
| isoCirc | 0.528 | 0.576 | 0.403 |
| circNICK-lrs | 0.739 | 0.767 | 0.836 |
| circFL-seq | 0.765 | 0.785 | 0.693 |
| `high_confidence` | 0.972 | 0.969 | 0.811 |
| `balanced_precision` | 0.973 | 0.966 | 0.828 |
| `balanced_recall` | 0.969 | 0.961 | 0.865 |
| `discovery` | 0.970 | 0.962 | 0.855 |

Every nanocirc tier correlates with ground truth expression more strongly than any single raw tool, on every protocol (r ≥ 0.81 vs raw tools' r 0.40-0.94). isoCirc has reached the weakest quantification of any tool/tier on every protocol (r 0.40-0.58), despite it's high precision on circRNA detection. CIRI-long is strongest quantifier among tools, but it fails to quantify well on isoCirc protocol.

<img src="images/benchmark/equal_sensitivity_by_expression_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_circfl.png" width="100%"/>
<img src="images/benchmark/equal_sensitivity_by_expression_isocirc.png" width="100%"/>

Detection sensitivity by ground-truth expression tertile (low/mid/high read count):

| Tier | Protocol | Low | Mid | High |
| --- | --- | --- | --- | --- |
| `discovery` | CIRI-long | 32.5% | 38.3% | 48.0% |
| `discovery` | circFL-seq | 40.6% | 47.1% | 56.7% |
| `discovery` | isoCirc | 13.9% | 21.3% | 34.2% |
| `balanced_recall` | CIRI-long | 31.7% | 37.8% | 47.3% |
| `balanced_recall` | circFL-seq | 40.0% | 46.3% | 56.3% |
| `balanced_recall` | isoCirc | 13.0% | 20.2% | 33.2% |
| `balanced_precision` | CIRI-long | 14.2% | 21.2% | 32.2% |
| `balanced_precision` | circFL-seq | 18.3% | 26.0% | 39.1% |
| `balanced_precision` | isoCirc | 3.4% | 7.6% | 16.8% |
| `high_confidence` | CIRI-long | 2.7% | 6.0% | 13.6% |
| `high_confidence` | circFL-seq | 3.8% | 8.2% | 16.7% |
| `high_confidence` | isoCirc | 0.1% | 0.4% | 2.6% |

Sensitivity grows with expression at every tier and protocol (more reads means more chances to be called). `balanced_recall` is close to `discovery` (within 1-2 points at every tertile), while `balanced_precision` is in between `discovery` and `high_confidence`. circFL-seq protocol allows the highest sensitivity at every tertile and tier; isoCirc protocol the lowest, most severely at `high_confidence`, where even highly-expressed circRNAs are barely detected (2.6% vs circFL-seq's 16.7%).

<img src="images/benchmark/equal_reproducibility_pairs_ciri_long.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_circfl.png" width="100%"/>
<img src="images/benchmark/equal_reproducibility_pairs_isocirc.png" width="100%"/>

Pairwise read-count agreement across every tool/tier combination, ground truth included. nanocirc's own tiers agree with ground truth better, than for any other caller and agree with each other almost perfectly, since they share the same remap-based quantification regardless of which tool found the locus. Agreement of raw tool calls with ground truth is looser and protocol-dependent; isoCirc shows the widest scatter against every other caller and ground truth on every protocol, consistent with its weak expression correlation above, CIRI-long shows the best agreement with ground truth among single tools on 2 out of 3 protocols.

</details>

**Protocol results:**

- **circFL-seq protocol generates reads that work well with all tools**, at every tier: highest recall, highest F1, and the highest BSJ-confirmation rate of the 3 protocols (64.1%).
- **CIRI-long protocol reads are also compatible with all tools**, with a lower BSJ-confirmation rate (54.8%) and lower recall at every tier than circFL-seq.
- **isoCirc protocol reads are the hardest for the tools to work with**, despite carrying the exact same ~9,800 underlying circRNAs. Its recall is the lowest of the 3 protocols at every tier, but precision remains comparable or higher to other wet-lab approaches. CIRI-long is almost incompatible with this wet lab protocol, and that was verified by running CIRI-long on real isoCIRC dataset as well - the tool never managed to finish the analysis on a full dataset.

Because the circRNA catalog, its abundances, and the simulated read amount were identical across all 3 comparisons, this difference comes entirely from how each protocol's wet-lab induced error profile and read length characteristics change the reads NanoSim generates from them.

**Recommendations**: if your protocol choice is flexible, circFL-seq or CIRI-long library prep will give nanocirc's tiers meaningfully more usable recall than isoCirc protocol. If you already generated your reads and are working with isoCirc data, you can prioritise balanced_recall mode which reaches high precision on this read type and expect potential incopatibility with CIRI-long tool.

<details>
<summary><strong>Statistical significance</strong></summary>

### Statistical significance

We tested whether the protocol-to-protocol differences above are real or could be noise. Each pair matches every caller/circRNA-type/run point between the two protocols, using these statistical methods (used the same way for every "statistical significance" comparison in this document):

- **Linear regression (slope, intercept, R²)**: fits one side's value from the other's (value = intercept + slope × other value). Slope near 1 with a small intercept means the two sides are close to equal plus a fixed shift. R² is how well that line predicts the value, whatever the size of the shift.
- **Kendall's tau-b**: ranks the 8 callers best to worst on each side, then checks how well the two rankings agree.
- **Mantel-Haenszel odds ratio + Breslow-Day test**: uses the raw TP/FP/FN counts directly, stratified by circRNA type. The odds ratio is how many times more likely a call is to be correct, or a ground-truth circRNA is recovered, on one side versus the other. Breslow-Day then checks if that odds ratio is the same size for every type; p < 0.05 means the effect size depends on circRNA type.
- **Wilcoxon signed-rank**: tests if one side's values are systematically higher than the other's, without assuming a normal distribution.
- **Cohen's d (paired)**: how big that gap is, in standard deviations. ~0.2 is small, ~0.5 medium, ~0.8+ large.

The two isoCirc-protocol pairs (CIRI-long vs isoCirc, circFL-seq vs isoCirc) are not included yet: this dataset's isoCirc-protocol run is missing CIRI-long output due to a pipeline default-drift regression (a rerun with CIRI-long explicitly re-enabled is in progress), so those two comparisons will be added once it completes.

#### CIRI-long vs circFL-seq

##### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | CIRI-long | circFL-seq | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 114 | 0.720 | 0.731 | +0.011 | 0.993 | 0.016 | 0.987 | 0.009 | 0.353 |
| Recall | 120 | 0.221 | 0.275 | +0.054 | 1.104 | 0.031 | 0.974 | 2×10⁻¹⁹ | 1.230 |
| F1 | 120 | 0.282 | 0.337 | +0.054 | 1.071 | 0.034 | 0.965 | 2×10⁻¹⁹ | 1.183 |

Precision fits an almost perfect 1:1 line (R² = 0.99, slope ≈ 1) but the gap itself is small (+0.011) and only modestly significant (p = 0.009, d = 0.35). Recall and F1 are both significantly and substantially higher on circFL-seq (p < 10⁻¹⁸, d ≈ 1.2), with a similarly tight fit (R² ≥ 0.96).

##### Per type, by metric (n=24 paired points per row pooling all 8 callers × 3 runs, except ciRNA precision n=21 and intergenic precision n=21: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | CIRI-long | circFL-seq | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 24 | 0.630 | 0.628 | −0.002 | 0.988 | 0.006 | 0.997 | 0.643 | −0.144 |
| eciRNA | Recall | 24 | 0.536 | 0.600 | +0.064 | 1.044 | 0.041 | 0.980 | 1×10⁻⁷ | 2.191 |
| eciRNA | F1 | 24 | 0.547 | 0.583 | +0.036 | 1.004 | 0.034 | 0.970 | 1×10⁻⁷ | 1.371 |
| EIciRNA | Precision | 24 | 0.635 | 0.655 | +0.020 | 0.979 | 0.033 | 0.970 | 0.060 | 0.485 |
| EIciRNA | Recall | 24 | 0.118 | 0.194 | +0.075 | 1.606 | 0.004 | 0.990 | 1×10⁻⁷ | 1.692 |
| EIciRNA | F1 | 24 | 0.188 | 0.281 | +0.093 | 1.401 | 0.018 | 0.985 | 1×10⁻⁷ | 2.049 |
| ciRNA | Precision | 21 | 0.879 | 0.905 | +0.026 | 0.737 | 0.257 | 0.971 | 0.003 | 0.828 |
| ciRNA | Recall | 24 | 0.145 | 0.208 | +0.064 | 1.451 | −0.001 | 0.999 | 6×10⁻⁵ | 1.120 |
| ciRNA | F1 | 24 | 0.222 | 0.296 | +0.074 | 1.306 | 0.006 | 0.999 | 6×10⁻⁵ | 1.364 |
| antisense | Precision | 24 | 0.693 | 0.683 | −0.010 | 0.985 | 0.001 | 0.999 | 0.003 | −0.726 |
| antisense | Recall | 24 | 0.156 | 0.183 | +0.028 | 1.166 | 0.002 | 0.990 | 2×10⁻⁴ | 1.009 |
| antisense | F1 | 24 | 0.241 | 0.272 | +0.031 | 1.116 | 0.003 | 0.991 | 2×10⁻⁴ | 0.977 |
| intergenic | Precision | 21 | 0.792 | 0.815 | +0.023 | 1.024 | 0.004 | 0.967 | 0.001 | 0.831 |
| intergenic | Recall | 24 | 0.151 | 0.189 | +0.038 | 1.246 | 0.001 | 1.000 | 6×10⁻⁵ | 1.024 |
| intergenic | F1 | 24 | 0.214 | 0.252 | +0.038 | 1.163 | 0.004 | 0.999 | 6×10⁻⁵ | 1.178 |

Every type's recall and F1 fit almost perfectly (R² ≥ 0.97) and circFL-seq is significantly ahead on both, for every type (p ≤ 2×10⁻⁴). Precision is smaller in magnitude and mixed in direction: eciRNA and antisense slightly favor CIRI-long, while ciRNA and intergenic significantly favor circFL-seq (p ≤ 0.003); EIciRNA shows no significant precision difference (p = 0.06).

##### Does caller ranking transfer?

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 8 | 0.929 (p = 4×10⁻⁴) | 0.976 (p = 3×10⁻⁵) |
| Recall | 8 | 1.000 (p = 5×10⁻⁵) | 1.000 (p = 0) |
| F1 | 8 | 1.000 (p = 5×10⁻⁵) | 1.000 (p = 0) |

Ranking transfers essentially perfectly. By F1:

| Caller | F1 CIRI-long | F1 circFL-seq | Rank CIRI-long | Rank circFL-seq |
| --- | --- | --- | --- | --- |
| balanced_recall | 0.470 | 0.558 | 1 | 1 |
| discovery | 0.457 | 0.538 | 2 | 2 |
| CIRI-long | 0.375 | 0.447 | 3 | 3 |
| balanced_precision | 0.304 | 0.364 | 4 | 4 |
| circFL-seq | 0.287 | 0.334 | 5 | 5 |
| isoCirc | 0.207 | 0.261 | 6 | 6 |
| high_confidence | 0.108 | 0.136 | 7 | 7 |
| circNICK-lrs | 0.051 | 0.055 | 8 | 8 |

All eight callers have on the exact same F1 rank on both protocols.

##### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=9)

| Tool | Metric | n | CIRI-long | circFL-seq | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.694 | 0.731 | +0.038 | 0.914 | 0.098 | 0.925 | 1×10⁻⁴ | 1.218 |
| CIRI-long | Recall | 15 | 0.293 | 0.360 | +0.067 | 0.954 | 0.081 | 0.947 | 6×10⁻⁵ | 1.563 |
| CIRI-long | F1 | 15 | 0.375 | 0.447 | +0.072 | 0.801 | 0.147 | 0.807 | 6×10⁻⁵ | 1.326 |
| isoCirc | Precision | 15 | 0.891 | 0.885 | −0.006 | 1.083 | −0.080 | 0.943 | 0.272 | −0.214 |
| isoCirc | Recall | 15 | 0.136 | 0.181 | +0.045 | 1.229 | 0.014 | 0.996 | 6×10⁻⁵ | 1.271 |
| isoCirc | F1 | 15 | 0.207 | 0.261 | +0.054 | 1.113 | 0.031 | 0.992 | 6×10⁻⁵ | 1.799 |
| circNICK-lrs | Precision | 9 | 0.096 | 0.091 | −0.004 | 0.994 | −0.004 | 0.986 | 0.249 | −0.473 |
| circNICK-lrs | Recall | 15 | 0.057 | 0.065 | +0.008 | 1.137 | 0.000 | 1.000 | 0.028 | 0.563 |
| circNICK-lrs | F1 | 15 | 0.051 | 0.055 | +0.003 | 1.061 | 0.000 | 0.999 | 0.028 | 0.602 |
| circFL-seq | Precision | 15 | 0.643 | 0.634 | −0.009 | 0.983 | 0.002 | 0.997 | 0.107 | −0.483 |
| circFL-seq | Recall | 15 | 0.203 | 0.241 | +0.038 | 0.969 | 0.045 | 0.974 | 6×10⁻⁵ | 1.281 |
| circFL-seq | F1 | 15 | 0.287 | 0.334 | +0.048 | 0.974 | 0.055 | 0.963 | 6×10⁻⁵ | 1.150 |
| high_confidence | Precision | 15 | 0.931 | 0.936 | +0.005 | 0.923 | 0.077 | 1.000 | 0.024 | 0.786 |
| high_confidence | Recall | 15 | 0.072 | 0.093 | +0.021 | 1.179 | 0.008 | 0.996 | 0.002 | 0.826 |
| high_confidence | F1 | 15 | 0.108 | 0.136 | +0.028 | 1.117 | 0.015 | 0.992 | 0.002 | 0.977 |
| balanced_precision | Precision | 15 | 0.898 | 0.907 | +0.009 | 1.118 | −0.097 | 0.870 | 0.762 | 0.279 |
| balanced_precision | Recall | 15 | 0.218 | 0.271 | +0.053 | 1.087 | 0.034 | 0.992 | 6×10⁻⁵ | 1.933 |
| balanced_precision | F1 | 15 | 0.304 | 0.364 | +0.060 | 0.980 | 0.066 | 0.977 | 6×10⁻⁵ | 1.914 |
| balanced_recall | Precision | 15 | 0.699 | 0.725 | +0.026 | 0.853 | 0.129 | 0.916 | 0.041 | 0.776 |
| balanced_recall | Recall | 15 | 0.391 | 0.491 | +0.100 | 0.956 | 0.117 | 0.972 | 6×10⁻⁵ | 3.017 |
| balanced_recall | F1 | 15 | 0.470 | 0.558 | +0.088 | 0.718 | 0.221 | 0.896 | 6×10⁻⁵ | 1.956 |
| discovery | Precision | 15 | 0.660 | 0.680 | +0.021 | 0.947 | 0.056 | 0.953 | 0.095 | 0.603 |
| discovery | Recall | 15 | 0.398 | 0.497 | +0.099 | 0.936 | 0.124 | 0.971 | 6×10⁻⁵ | 2.820 |
| discovery | F1 | 15 | 0.457 | 0.538 | +0.082 | 0.619 | 0.256 | 0.783 | 6×10⁻⁵ | 1.638 |

Every caller's recall and F1 is significantly higher on circFL-seq (p ≤ 6×10⁻⁵) with a good fit (R² ≥ 0.78). For precision CIRI-long, high_confidence, ciRNA-heavy balanced_precision, and balanced_recall show a small but significant circFL-seq precision edge (p ≤ 0.041), while isoCirc, circNICK-lrs, circFL-seq-seq, and discovery show no significant precision difference at all between the two protocols.

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 1.15 | 1.09–1.21 | 4×10⁻⁸ | 8×10⁻⁶ |
| CIRI-long | Recall | 5 | 1.43 | 1.38–1.49 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 0.88 | 0.79–0.98 | 0.018 | 0.075 |
| isoCirc | Recall | 5 | 1.52 | 1.45–1.60 | <10⁻³⁰⁰ | 5×10⁻⁴ |
| circNICK-lrs | Precision | 5 | 1.00 | 0.93–1.07 | 0.910 | 0.096 |
| circNICK-lrs | Recall | 5 | 1.19 | 1.10–1.28 | 6×10⁻⁶ | 0.794 |
| circFL-seq | Precision | 5 | 0.93 | 0.87–1.01 | 0.084 | 0.988 |
| circFL-seq | Recall | 5 | 1.31 | 1.26–1.37 | <10⁻³⁰⁰ | 4×10⁻¹⁵ |
| high_confidence | Precision | 5 | 1.12 | 0.96–1.30 | 0.147 | 0.897 |
| high_confidence | Recall | 5 | 1.44 | 1.35–1.54 | <10⁻³⁰⁰ | 1×10⁻⁵ |
| balanced_precision | Precision | 5 | 0.92 | 0.85–1.00 | 0.058 | 2×10⁻⁴ |
| balanced_precision | Recall | 5 | 1.45 | 1.39–1.52 | <10⁻³⁰⁰ | 3×10⁻⁴ |
| balanced_recall | Precision | 5 | 1.11 | 1.06–1.16 | 5×10⁻⁶ | 4×10⁻⁸ |
| balanced_recall | Recall | 5 | 1.61 | 1.55–1.67 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 1.05 | 1.01–1.10 | 0.010 | 2×10⁻¹⁰ |
| discovery | Recall | 5 | 1.60 | 1.55–1.66 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

Every caller has significantly higher odds of recovering a ground-truth circRNA (recall) on circFL-seq (OR 1.2–1.6). Precision odds are mixed and mostly weaker: isoCirc actually favors CIRI-long (OR = 0.88, p = 0.018), circNICK-lrs/circFL-seq/high_confidence show no significant precision difference, and the rest favor circFL-seq modestly (OR 1.05–1.15). 5 of 16 caller/outcome combinations have a protocol effect that's uniform across circRNA types (circNICK-lrs's precision and recall, isoCirc's precision, high_confidence's precision, circFL-seq's precision), the most homogeneity of any comparison in this document; the other 11 combinations' boost size still depends on which circRNA type it is.

**Conclusion.** Correlation is excellent for both metrics (aggregate R² ≥ 0.97, essentially every type and caller individually significant on recall/F1). Caller ranking transfers perfectly (Kendall's tau-b = 1.00 on recall/F1, all 8 callers same rank by F1). Precision differences are real but small and caller/type-dependent, with several callers and types showing no significant precision difference between protocols while circFL-seq is better for recall/F1.

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

CIRI-long is not run by default on this protocol (see below), so it has no standalone FP characterization here.

| FP read is... | CIRI-long | circFL-seq | isoCirc | circNICK-lrs |
| --- | --- | --- | --- | --- |
| Not from a circ read at all (linear read miscalled) | not run | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.2 ± 0.1 |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | not run | 28.9 ± 0.2 | 91.8 ± 0.9 | 1.9 ± 0.1 |
| Real circ read, correct strand, pure boundary noise | not run | 60.2 ± 0.3 | 0.0 ± 0.0 | 18.5 ± 1.2 |
| Real circ read, but circRNA reported at a different/distant locus | not run | 10.9 ± 0.3 | 8.2 ± 0.9 | 79.4 ± 1.3 |
| FP locus overlaps a genomic repeat | not run | 52.4 ± 0.8 | 22.4 ± 1.4 | 57.4 ± 0.5 |
| n FP loci (pooled, 3 runs) | not run | 2,269 | 354 | 10,995 |
| n supporting reads (pooled, 3 runs) | not run | 13,832 | 1,103 | 29,882 |

Many FP come from `antisense`-type circRNA, where the simulator deliberately sets the circRNA's strand opposite to its host gene; since strand agreement is required by every matching strategy, no boundary-tolerance setting can rescue these. That pattern is present for CIRI-long and isoCirc on all 3 protocols (77-94% wrong-strand), and for circFL-seq on the CIRI-long and circFL-seq protocols (73-77%), but circFL-seq's FP composition changes on isoCirc-protocol reads: wrong-strand drops to 28.9% and "correct strand, pure boundary noise" goes up to 60.2%, the opposite of its own pattern everywhere else and of every other tool on this protocol. circNICK-lrs's FPs are structurally different from the other 3 tools throughout: mostly boundary noise and distant-locus calls (66-79%), rarely a strand flip (2-3%), and consistently the most repeat-associated (57-63% vs 18-38% for the other tools).

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
| Share of GT catalog | 77.1% | 22.9% |
| Mean mature length (bp) | 2,486 | 583 |
| Mean rolling-circle copy number | 7.6 | 10.0 |
| Median simulated read count | 49 | 94 |
| Canonical GT/AG splice site at BSJ | 63.7% | 84.8% |
| No gene overlap at all (intergenic) | 21.4% | 10.6% |
| Falls in a gene with ≥5 isoforms | 72.4% | 81.6% |
| Overlaps a genomic repeat | 49.8% | 56.7% |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 43.4% | 56.6% |
| `antisense` | 90.7% | 9.3% |
| `ciRNA` (purely intronic) | 85.2% | 14.8% |
| `EIciRNA` (exon+intron) | 79.2% | 20.8% |
| `intergenic` | 87.9% | 12.1% |

circFL-seq protocol has the largest recovery of the 3 (48.0% of the GT catalog found by at least one caller), CIRI-long protocol is close behind (39.4%), and isoCirc protocol has by far the smallest (22.9%, barely more than half of CIRI-long protocol's). The missed-vs-recovered pattern is consistent across all 3 protocols regardless of that gap: missed circRNAs are longer (~5-9x), lower-expressed, and use non-canonical splice sites more often than recovered ones. `eciRNA` is always the easiest biotype to recover in all wet-lab protocols. On isoCirc protocol every biotype's recall collapses further than on the other 2 protocols, most severely for `antisense` (9.3% recovered, versus 33.7% on CIRI-long protocol and 38.6% on circFL-seq protocol), also due to the same wrong-strand FP pattern documented above rather than due to a different reason. Repeat overlap is not a meaningful driver of recovery on any protocol.

</details>

<details>
<summary><strong>isoCirc protocol: with CIRI-long off by default</strong></summary>

#### isoCirc protocol: with CIRI-long off by default

`--wet_lab isocirc` defaults `--run_cirilong` to `false` (see [Which tools work with my wet-lab protocol and species?](#which-tools-work-with-my-wet-lab-protocol-and-species) above). The isoCirc protocol column above used all 4 tools to estimate real performance of each tool; but we show here the results on the same dataset with CIRI-long left out, next to the all-4-tools numbers for comparison. 

<img src="images/benchmark/pr_scatter_equal_isocirc_no_cirilong.png" width="100%"/>

Mean precision (P) / recall (R) / F1 across the same 3 runs, exon-based matching, plus each individual tool's own raw calls:

| Tier / tool | CIRI-long on (old) | CIRI-long off (new default) |
| --- | --- | --- |
| isoCirc alone | P 0.87 / R 0.08 / F1 0.15 | P 0.87 / R 0.08 / F1 0.15 |
| circFL-seq alone | P 0.74 / R 0.23 / F1 0.36 | P 0.74 / R 0.23 / F1 0.36 |
| CIRI-long alone | P 0.77 / R 0.03 / F1 0.05 | not run |
| circNICK-lrs alone | P 0.22 / R 0.10 / F1 0.14 | P 0.22 / R 0.10 / F1 0.14 |
| `discovery` | P 0.54 / R 0.28 / F1 0.37 | P 0.66 / R 0.26 / F1 0.37 |
| `balanced_recall` | P 0.81 / R 0.25 / F1 0.38 | P 0.80 / R 0.24 / F1 0.38 |
| `balanced_precision` | P 0.85 / R 0.09 / F1 0.17 | P 0.85 / R 0.09 / F1 0.17 |
| `high_confidence` | P 0.92 / R 0.02 / F1 0.04 | P 0.87 / R 0.01 / F1 0.02 |

CIRI-long was the weakest tool on this protocol (F1 0.05, far behind its performance on the CIRI-long/circFL-seq protocols).

Dropping CIRI-long costs essentially nothing at any of the 4 tiers in this benchmark: P/R/F1 are identical or within 0.01 of each other with CIRI-long on vs off. Against the real-world risk documented above (CIRI-long hanging for 24+ days on real isoCirc-protocol data), this is why CIRI-long off is now the default for this protocol. If you wish to try to run it on isoCirc data despite the risk, you can turn it on yourself with `--run_cirilong true`.

circNICK-lrs wet-lab protocol reads were not tested, as all 3 other tools expect rolling-circle generated data with multiple copies of circRNA per read, making circNICK-lrs protocol compatible only with circNICK-lrs tool, thus it wouldn't benefit from nanocirc's consensus voting algorithm and it defeats the purpose of the pipeline. You still can run circNICK-lrs inside the nanocirc pipeline and benefit from nanocirc's upstream patches to circNICK-lrs tool and nanocirc quantification approach, but without confidence scoring or a sensible Precision-Recall curve.

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
| `discovery` | P 0.67 / R 0.37 / F1 0.48 | P 0.74 / R 0.54 / F1 0.62 |
| `balanced_recall` | P 0.74 / R 0.36 / F1 0.48 | P 0.85 / R 0.52 / F1 0.65 |
| `balanced_precision` | P 0.85 / R 0.20 / F1 0.33 | P 0.94 / R 0.36 / F1 0.52 |
| `high_confidence` | P 0.89 / R 0.07 / F1 0.13 | P 0.98 / R 0.13 / F1 0.22 |

</details>

<details>
<summary><strong>Ground truth coverage</strong></summary>

### Ground truth coverage

<img src="images/benchmark/gt_coverage_stacked_human.png" width="100%"/>
<img src="images/benchmark/gt_coverage_stacked_mouse.png" width="100%"/>

Mouse has higher ground truth recovery than human at every tier (e.g. `discovery` R 0.54 vs human's 0.37, from the table above), the same species gap the summary table shows.

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
| Human | 73.2% | 17.2% | 25.8% | 28.6% | 20.8% |
| Mouse | 88.8% | 38.9% | 53.5% | 54.0% | 15.9% |

Mouse recovers every type better than human except intergenic, where human is slightly ahead (20.8% vs 15.9%), the one type/species combination where human derived circRNAs win. eciRNA is the easiest type on both species; intergenic is hardest on mouse, EIciRNA hardest on human.

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

| Tier | Human | Mouse |
| --- | --- | --- |
| `discovery` | 31.0% | 41.9% |
| `balanced_recall` | 29.4% | 40.6% |
| `balanced_precision` | 20.4% | 27.5% |
| `high_confidence` | 6.9% | 10.1% |

Mouse resolves more isoforms than human at every tier, consistent with mouse's higher overall recall in this benchmark. Both species are far below the wet-lab-protocol benchmark's exact-isoform rates above (e.g. `discovery` 31.0%/41.9% here vs 61.8-67.6% on CIRI-long/circFL-seq protocol), since this dataset uses a more isoform-family-rich catalog (422/912 isoforms tested here vs 102 there).

<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_mouse.png" width="100%"/>

<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_90_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_90_mouse.png" width="100%"/>

Per family plot instead of per isoform: did the caller get every isoform in the family right, some of them, or none. % of isoform families fully resolved, at 99% reciprocal overlap:

| Tier | Human | Mouse |
| --- | --- | --- |
| `discovery` | 16.3% | 26.0% |
| `balanced_recall` | 14.4% | 23.5% |
| `balanced_precision` | 7.2% | 10.8% |
| `high_confidence` | 1.4% | 4.7% |

Getting one isoform right while missing or misclassifying others is the most common outcome; fully resolving an isoform family is rare for every tier on both species, and `high_confidence` barely ever manages it (1.4%/4.7%). Mouse is consistently ahead of human here too.

</details>

<details>
<summary><strong>Expression quantification and reproducibility</strong></summary>

### Expression quantification and reproducibility

<img src="images/benchmark/expression_correlation_human.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse.png" width="100%"/>

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
| circFL-seq | 0.787 | 0.805 |
| `high_confidence` | 0.967 | 0.930 |
| `balanced_precision` | 0.974 | 0.947 |
| `balanced_recall` | 0.967 | 0.936 |
| `discovery` | 0.967 | 0.936 |

isoCirc has the weakest quantification of any tool/tier on both species (r 0.55-0.67), same as on every wet-lab protocol above. Every nanocirc tier correlates more strongly with ground truth than any raw tool on both species (r ≥ 0.93). Unlike detection, where mouse leads, quantification quality is close between species and `balanced_precision` is actually slightly stronger on human (0.974 vs 0.947).

#### Isoform-level expression quantification

The correlation above evaluates BSJ level expression, so it only shows if a locus's read count is right, not per-isoform count. The three plots below show the same correlation, but matched at the `full`/`strict`/`relaxed` isoform tiers described in [Isoform detection](#isoform-detection-1) above. This shows how well tools estimate isoform level expression, not just locus level. Pearson r at the strict (99% overlap) isoform tier; n is small since only isoform-family loci qualify:

| Caller | Human r (n) | Mouse r (n) |
| --- | --- | --- |
| CIRI-long | 0.890 (77) | 0.818 (281) |
| isoCirc | 0.580 (54) | 0.515 (189) |
| circNICK-lrs | not enough calls (13) | 0.660 (51) |
| circFL-seq | 0.518 (75) | 0.665 (212) |
| `high_confidence` | 0.877 (44) | 0.845 (188) |
| `balanced_precision` | 0.901 (90) | 0.874 (304) |
| `balanced_recall` | 0.923 (122) | 0.870 (378) |
| `discovery` | 0.860 (125) | 0.870 (379) |

`balanced_recall` has the strongest isoform-level correlation on human (r=0.923); `balanced_precision` is the best on mouse (r=0.874), but the four nanocirc tiers are close to each other on both species. circFL-seq and isoCirc are the weakest quantifiers at this level too, same pattern as locus-level expression above.

<img src="images/benchmark/expression_correlation_human_full.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_full.png" width="100%"/>
<img src="images/benchmark/expression_correlation_human_strict.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_strict.png" width="100%"/>
<img src="images/benchmark/expression_correlation_human_relaxed.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse_relaxed.png" width="100%"/>

<img src="images/benchmark/sensitivity_by_expression_human.png" width="100%"/>
<img src="images/benchmark/sensitivity_by_expression_mouse.png" width="100%"/>

Detection sensitivity broken down by ground-truth expression tertile: low-expressed circRNAs are structurally harder to detect correctly than high-expressed ones, for every tool.

| Tier | Species | Low | Mid | High |
| --- | --- | --- | --- | --- |
| `discovery` | Human | 15.7% | 37.3% | 53.5% |
| `discovery` | Mouse | 26.7% | 55.3% | 68.6% |
| `high_confidence` | Human | 0.8% | 5.8% | 15.3% |
| `high_confidence` | Mouse | 2.4% | 11.4% | 22.8% |

Mouse is more sensitive than human at every tertile and tier shown, widest at `discovery`'s low tertile (26.7% vs 15.7%, a 70% relative gap) and narrowest at `high_confidence`'s low tertile, where both species are already close to floor.

<img src="images/benchmark/reproducibility_pairs_human.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse.png" width="100%"/>

Pairwise read-count agreement across every tool/tier combination, ground truth included. As with the wet-lab-protocol benchmark above, nanocirc's own tiers agree with each other almost perfectly (diagonal points in every tier-vs-tier panel) since they share the same remap-based quantification; agreement against raw tool calls is looser, isoCirc again showing the widest scatter, consistent with its weaker Pearson r above.

Same reproducibility matrix, matched at the `full`/`strict`/`relaxed` isoform tiers instead of structural comparison, so agreement is evaluated on isoform level.

<img src="images/benchmark/reproducibility_pairs_human_full.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_full.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_human_strict.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_strict.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_human_relaxed.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse_relaxed.png" width="100%"/>

</details>

**What we see:**

- **Mouse outperformed human at every tier**, on both precision and recall. The difference is largest at `balanced_precision` tier (F1 0.50 vs 0.32) and smallest at `discovery` tier (F1 0.59 vs 0.47).

- **Quantification differences between species are not as big**, on both precision and recall. Quantification performance is close between species and `balanced_precision` is actually slightly stronger on human (0.974 vs 0.947).

**Recommendations**: species-level difference for detection are present, especially for recall and F1. Published tier-level precision/recall for human or mouse are not as guaranteed to transfer across other species. Re-benchmark on your own species/catalog using nanocirc-simulator and nanocirc-benchmark and turn off circNICK-lrs for any species but human or mouse (since it's incompatible with any other species).

<details>
<summary><strong>Statistical significance</strong></summary>

### Statistical significance

Same statistical methods as the [wet-lab-protocol comparison](#statistical-significance) above (linear regression, Kendall's tau-b, Mantel-Haenszel + Breslow-Day, Wilcoxon signed-rank, Cohen's d), matching every caller/circRNA-type/run point between human and mouse.

#### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | Human | Mouse | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 108 | 0.725 | 0.734 | +0.009 | 0.935 | 0.057 | 0.836 | 0.436 | 0.092 |
| Recall | 120 | 0.187 | 0.290 | +0.103 | 1.213 | 0.063 | 0.852 | 2×10⁻¹⁶ | 0.977 |
| F1 | 120 | 0.248 | 0.361 | +0.113 | 1.240 | 0.054 | 0.862 | 6×10⁻¹⁶ | 1.023 |

Recall and F1 are well predicted by human's value (R² 0.85–0.86) with a slope above 1: mouse's advantage is higher for callers that are stronger, not frowing uniformly for all callers. Precision fits well (R² = 0.84, slope near 1) but shows no systematic species gap at all (Wilcoxon p = 0.44): mouse is not reliably more precise than human once every caller and type is pooled together.

#### Per type, by metric (n=24 paired points per row pooling all 8 callers × 3 runs, except ciRNA precision n=21, antisense precision n=18, and intergenic precision n=21: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | Human | Mouse | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 24 | 0.651 | 0.748 | +0.096 | 1.046 | 0.066 | 0.945 | 1×10⁻⁷ | 1.699 |
| eciRNA | Recall | 24 | 0.495 | 0.619 | +0.124 | 1.111 | 0.069 | 0.935 | 1×10⁻⁷ | 2.196 |
| eciRNA | F1 | 24 | 0.532 | 0.657 | +0.125 | 1.206 | 0.016 | 0.938 | 1×10⁻⁷ | 2.280 |
| EIciRNA | Precision | 24 | 0.629 | 0.652 | +0.023 | 0.974 | 0.039 | 0.827 | 0.747 | 0.226 |
| EIciRNA | Recall | 24 | 0.103 | 0.230 | +0.127 | 2.287 | −0.006 | 0.961 | 1×10⁻⁷ | 1.725 |
| EIciRNA | F1 | 24 | 0.168 | 0.321 | +0.153 | 1.855 | 0.009 | 0.924 | 2×10⁻⁶ | 1.895 |
| ciRNA | Precision | 21 | 0.862 | 0.925 | +0.063 | 0.562 | 0.441 | 0.984 | 2×10⁻⁴ | 1.331 |
| ciRNA | Recall | 24 | 0.114 | 0.242 | +0.128 | 2.030 | 0.011 | 0.995 | 6×10⁻⁵ | 1.219 |
| ciRNA | F1 | 24 | 0.182 | 0.336 | +0.153 | 1.658 | 0.033 | 0.990 | 6×10⁻⁵ | 1.498 |
| antisense | Precision | 18 | 0.774 | 0.720 | −0.054 | 0.939 | −0.006 | 0.997 | 4×10⁻⁴ | −1.860 |
| antisense | Recall | 24 | 0.135 | 0.278 | +0.143 | 1.984 | 0.011 | 0.863 | 2×10⁻⁴ | 0.976 |
| antisense | F1 | 24 | 0.216 | 0.356 | +0.141 | 1.603 | 0.010 | 0.894 | 2×10⁻⁴ | 0.952 |
| intergenic | Precision | 21 | 0.742 | 0.635 | −0.107 | 0.896 | −0.029 | 0.853 | 1×10⁻⁵ | −1.691 |
| intergenic | Recall | 24 | 0.089 | 0.083 | −0.006 | 0.714 | 0.020 | 0.922 | 0.434 | −0.191 |
| intergenic | F1 | 24 | 0.139 | 0.133 | −0.006 | 0.718 | 0.033 | 0.905 | 0.794 | −0.138 |

Mouse is significantly ahead on recall and F1 for eciRNA, EIciRNA, ciRNA, and antisense. intergenic is the exception: mouse is not significantly different from human on either recall or F1 for this type (p = 0.43 and 0.79). For precision species relationship flips by type: eciRNA and ciRNA favor mouse (p ≤ 2×10⁻⁴), EIciRNA shows no significant difference (p = 0.75), and antisense and intergenic actually favor human (negative differences, p ≤ 4×10⁻⁴, the only rows in this whole table where human is significantly ahead of mouse).

#### Does caller ranking transfer?

We averaged each caller's human and mouse values across type and run first, ranked the 8 callers by that average on each side, then compared the two rankings with Kendall's tau-b:

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 8 | 0.929 (p = 4×10⁻⁴) | 0.976 (p = 3×10⁻⁵) |
| Recall | 8 | 0.929 (p = 4×10⁻⁴) | 0.976 (p = 3×10⁻⁵) |
| F1 | 8 | 0.929 (p = 4×10⁻⁴) | 0.976 (p = 3×10⁻⁵) |

Ranking transfers significantly on every metric. By F1:

| Caller | F1 human | F1 mouse | Rank human | Rank mouse |
| --- | --- | --- | --- | --- |
| balanced_recall | 0.406 | 0.572 | 1 | 1 |
| discovery | 0.396 | 0.545 | 2 | 2 |
| CIRI-long | 0.319 | 0.442 | 3 | 3 |
| balanced_precision | 0.263 | 0.440 | 4 | 4 |
| circFL-seq | 0.251 | 0.327 | 5 | 6 |
| isoCirc | 0.188 | 0.343 | 6 | 5 |
| high_confidence | 0.102 | 0.160 | 7 | 7 |
| circNICK-lrs | 0.056 | 0.055 | 8 | 8 |

Six of the eight callers reach the exact same F1 rank on both species. Only circFL-seq and isoCirc swap places (circFL-seq ranks 5th on human but 6th on mouse, isoCirc the reverse).

#### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=6 and high_confidence precision n=12)

| Tool | Metric | n | Human | Mouse | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.687 | 0.791 | +0.104 | 0.389 | 0.524 | 0.350 | 0.003 | 1.078 |
| CIRI-long | Recall | 15 | 0.240 | 0.338 | +0.097 | 1.053 | 0.085 | 0.833 | 6×10⁻⁴ | 1.124 |
| CIRI-long | F1 | 15 | 0.319 | 0.442 | +0.123 | 1.132 | 0.081 | 0.713 | 9×10⁻⁴ | 1.153 |
| isoCirc | Precision | 15 | 0.877 | 0.844 | −0.032 | 0.880 | 0.073 | 0.478 | 0.135 | −0.426 |
| isoCirc | Recall | 15 | 0.121 | 0.255 | +0.133 | 1.581 | 0.063 | 0.866 | 6×10⁻⁵ | 1.148 |
| isoCirc | F1 | 15 | 0.188 | 0.343 | +0.155 | 1.362 | 0.087 | 0.908 | 6×10⁻⁵ | 1.429 |
| circNICK-lrs | Precision | 6 | 0.167 | 0.131 | −0.036 | 4.369 | −0.598 | 0.906 | 0.438 | −0.709 |
| circNICK-lrs | Recall | 15 | 0.056 | 0.060 | +0.004 | 1.036 | 0.002 | 0.999 | 0.028 | 0.774 |
| circNICK-lrs | F1 | 15 | 0.056 | 0.055 | −0.000 | 1.037 | −0.002 | 0.994 | 0.600 | −0.047 |
| circFL-seq | Precision | 15 | 0.624 | 0.603 | −0.020 | 1.053 | −0.053 | 0.946 | 0.972 | −0.237 |
| circFL-seq | Recall | 15 | 0.174 | 0.236 | +0.062 | 1.081 | 0.048 | 0.918 | 0.009 | 1.055 |
| circFL-seq | F1 | 15 | 0.251 | 0.327 | +0.076 | 1.074 | 0.057 | 0.886 | 0.013 | 0.924 |
| high_confidence | Precision | 12 | 0.925 | 0.902 | −0.023 | 0.475 | 0.462 | 0.240 | 0.369 | −0.291 |
| high_confidence | Recall | 15 | 0.065 | 0.109 | +0.044 | 1.496 | 0.012 | 0.995 | 6×10⁻⁵ | 0.759 |
| high_confidence | F1 | 15 | 0.102 | 0.160 | +0.058 | 1.347 | 0.023 | 0.990 | 6×10⁻⁵ | 0.913 |
| balanced_precision | Precision | 15 | 0.888 | 0.872 | −0.016 | 0.322 | 0.586 | 0.088 | 0.421 | −0.197 |
| balanced_precision | Recall | 15 | 0.186 | 0.332 | +0.146 | 1.128 | 0.122 | 0.858 | 6×10⁻⁵ | 1.563 |
| balanced_precision | F1 | 15 | 0.263 | 0.440 | +0.177 | 1.074 | 0.157 | 0.870 | 6×10⁻⁵ | 1.961 |
| balanced_recall | Precision | 15 | 0.690 | 0.740 | +0.050 | 0.662 | 0.283 | 0.370 | 0.095 | 0.412 |
| balanced_recall | Recall | 15 | 0.323 | 0.492 | +0.168 | 0.999 | 0.169 | 0.740 | 9×10⁻⁴ | 1.382 |
| balanced_recall | F1 | 15 | 0.406 | 0.572 | +0.166 | 1.073 | 0.137 | 0.593 | 9×10⁻⁴ | 1.277 |
| discovery | Precision | 15 | 0.650 | 0.664 | +0.014 | 0.837 | 0.120 | 0.722 | 0.890 | 0.154 |
| discovery | Recall | 15 | 0.331 | 0.502 | +0.171 | 1.010 | 0.168 | 0.756 | 9×10⁻⁴ | 1.410 |
| discovery | F1 | 15 | 0.396 | 0.545 | +0.149 | 1.121 | 0.101 | 0.561 | 9×10⁻⁴ | 1.199 |

CIRI-long, isoCirc, high_confidence, balanced_precision, balanced_recall, and discovery are all significantly higher on mouse for recall and F1 (p ≤ 9×10⁻⁴). circFL-seq and circNICK-lrs are the exceptions: circFL-seq's recall/F1 gap is smaller (p ≤ 0.013), and circNICK-lrs shows no significant recall or F1 difference between species at all. Only CIRI-long shows a significant species gap on precision (p = 0.003, mouse higher); every other caller's precision is statistically indistinguishable between human and mouse.

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 2.17 | 2.03–2.31 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| CIRI-long | Recall | 5 | 1.79 | 1.72–1.87 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 1.03 | 0.92–1.16 | 0.599 | <10⁻³⁰⁰ |
| isoCirc | Recall | 5 | 3.11 | 2.96–3.27 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circNICK-lrs | Precision | 5 | 0.90 | 0.83–0.97 | 0.004 | <10⁻³⁰⁰ |
| circNICK-lrs | Recall | 5 | 1.08 | 1.00–1.16 | 0.054 | 0.172 |
| circFL-seq | Precision | 5 | 1.06 | 0.97–1.15 | 0.214 | <10⁻³⁰⁰ |
| circFL-seq | Recall | 5 | 1.64 | 1.56–1.72 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| high_confidence | Precision | 5 | 1.88 | 1.57–2.26 | 1×10⁻¹¹ | 5×10⁻⁵ |
| high_confidence | Recall | 5 | 2.10 | 1.96–2.25 | <10⁻³⁰⁰ | 5×10⁻⁵ |
| balanced_precision | Precision | 5 | 1.53 | 1.39–1.68 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced_precision | Recall | 5 | 2.73 | 2.61–2.86 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced_recall | Precision | 5 | 1.59 | 1.51–1.68 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced_recall | Recall | 5 | 2.45 | 2.36–2.55 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 1.26 | 1.21–1.32 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Recall | 5 | 2.51 | 2.41–2.61 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

Every caller has significantly higher odds of recovering a ground-truth circRNA (recall) on mouse than on human (OR 1.1–3.1), and most also have higher odds of a correct call (precision), though isoCirc and circFL-seq show no significant precision difference (p = 0.60 and 0.21) and circNICK-lrs's precision actually favors human (OR = 0.90, p = 0.004). Only 1 of the 16 caller/outcome combinations (circNICK-lrs's recall) has a species effect that's is stable across circRNA types (Breslow-Day p = 0.17); every other combination's effect size depends on circRNA type.

**Conclusion.** Correlation is good for recall and F1, but weaker and inconsistent for precision, where half of callers and types show no significant species difference at all, and antisense/intergenic precision actually favors human. Caller ranking transfers strongly (Kendall's tau-b = 0.93 on every metric); only circFL-seq and isoCirc swap adjacent ranks. The size of the species gap depends on circRNA types (1 of 16 Mantel-Haenszel/Breslow-Day combinations homogeneous). So the recall/F1 species gap exists, the caller ranking is mostly preserved between the species, and precision across species is mostly caller- and type-dependent.

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

Other benchmarks use synthetic circRNAs generated by nanocirc-simulator and processed into reads by NanoSim. This benchmark evaluates performance using 2,000 database-derived circRNAs per type instead, sampled randomly from the human circBase/circAtlas union catalog. It uses the same NanoSim read simulation, expression distribution, detection pipeline, and quantification as every other dataset here. It shows the recovery of already-annotated circRNAs versus in silico de novo ones. It is important to note that these databases can contain circRNAs derived from short-read sequencing computational predictions that were never verified in a wet lab, are biased toward what tools can already detect (since that is how they were catalogued in the first place), and that some tools use the circRNA database itself as direct input for detection. It is expected that tools will perform better on database circRNAs, and there is no guarantee that every circRNA in the database is biologically true — but even with these limitations, it is a useful comparator and the closest approximation to biological ground truth we have.

<img src="images/benchmark/pr_scatter_database_human.png" width="100%"/>

<details>
<summary><strong>Database benchmark results</strong></summary>

### Database benchmark results

<img src="images/benchmark/pr_scatter_isoform_database_human.png" width="100%"/>

Same precision/recall, but at strict isoform-level matching (BSJ ±10bp and ≥99% exon overlap), the same criteria used for every comparison in this section.

<img src="images/benchmark/gt_coverage_stacked_database_human.png" width="100%"/>

Green (recovered) plus grey (missed) sum to 100% of ground truth; red stacks separately on top as false calls, scaled the same way.

<img src="images/benchmark/type_classification_database_human_isoform_level.png" width="100%"/>

Same TP/FP vs GT breakdown, split by circRNA biotype (`eciRNA`, `EIciRNA`, `ciRNA`, `antisense`, `intergenic`) instead of one bar per caller.

<img src="images/benchmark/isoform_detection_reciprocal_overlap_99_database_human.png" width="100%"/>
<img src="images/benchmark/isoform_family_completeness_reciprocal_overlap_99_database_human.png" width="100%"/>

Isoform-level detection accuracy and per-family completeness within shared-BSJ multi-isoform families.

<img src="images/benchmark/expression_correlation_database_human.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_database_human.png" width="100%"/>

Read-count correlation against ground truth, and pairwise cross-tool/tier agreement, per caller.

</details>

### Does the synthetic benchmark conclusions transfer to the database-derived dataset?

Database circRNAs are easier to detect than synthetic ones by design. They are known, already-catalogued circRNAs, already selected for detectability, and used as input or a built-in dependency by tools during circRNA calling. The in silico dataset instead tries to evaluate performance on novel, unannotated molecules: it asks "if the read evidence for this molecule is there, but no annotation hints exist, can the tools still recover it?" The database benchmark asks a different question instead: "of circRNAs that are already known, annotated, and previously detected by different tools, how many can a tool recover when given read evidence and, where possible, the same input database?"

So the gap between the two methods is expected, but it's hard to say how much of it is due to the simulation logic, the database being biased by pre-selection for detectability, the presence of circRNA structures that were already reported by the same tools (which may still be real TPs or FPs, since we don't have wet-lab validation data for each database entry and depositing one doesn't require it), the presence of a database input file that directly hints the tools, actual biology, or other factors. What we can measure is whether tool ranking and per-type detectability can be transferred between the two. 

We matched precision/recall/F1 between the in silico and database benchmarks at strict isoform-level matching (BSJ ±10bp and ≥99% exon overlap), separately for every caller, circRNA type, and run: up to 120 paired points (8 callers × 5 types × 3 runs; precision has 111, since a caller that made zero calls of a type in a run has undefined precision, not 0, so that cell is dropped). Computed by [`scripts/statistical_analysis.py`](https://gitlab.com/bingo-igdr/nano-circ), reusable for any two datasets scored with our benchmark (not just in silico vs database, also wet-lab protocol vs protocol and species vs species, both used elsewhere in this document).

Same statistical methods as the [wet-lab-protocol comparison](#statistical-significance) (linear regression, Kendall's tau-b, Mantel-Haenszel + Breslow-Day, Wilcoxon signed-rank, Cohen's d), applied here to synthetic vs database instead of protocol vs protocol.

<details>
<summary><strong>Aggregate, by metric</strong></summary>

#### Aggregate, by metric (n paired points pooling every caller × type × run)

| Metric | n | Synthetic | Database | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Precision | 111 | 0.679 | 0.828 | +0.150 | 0.654 | 0.384 | 0.611 | 1×10⁻¹⁷ | 0.850 |
| Recall | 120 | 0.187 | 0.326 | +0.139 | 1.084 | 0.123 | 0.621 | 9×10⁻¹⁹ | 0.851 |
| F1 | 120 | 0.248 | 0.415 | +0.168 | 1.164 | 0.127 | 0.691 | 6×10⁻¹⁹ | 1.051 |

Recall and F1 on database circRNAs are well predicted by their in silico counterpart values (R² 0.62–0.69), close to a straight 1:1 relationship (slope 1.08–1.16). Precision is less predictable (R² = 0.61) and doesn't scale evenly: callers with weak synthetic-side precision gain the most on database circRNAs, while callers already near the ceiling gain little (slope 0.65). All three metrics are significantly higher on database circRNAs (Wilcoxon p < 10⁻¹⁶, Cohen's d 0.85–1.05).

</details>

<details>
<summary><strong>Per type, by metric</strong></summary>

#### Per type, by metric (n=24 paired points per row pooling all 8 callers × 3 runs, except antisense precision, n=21, and intergenic precision, n=18: some caller × run cells had zero calls of that type, so precision can't be defined and those cells are dropped)

| Type | Metric | n | Synthetic | Database | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| eciRNA | Precision | 24 | 0.651 | 0.804 | +0.153 | 0.835 | 0.260 | 0.858 | 1×10⁻⁷ | 1.831 |
| eciRNA | Recall | 24 | 0.495 | 0.550 | +0.055 | 1.091 | 0.010 | 0.987 | 1×10⁻⁷ | 1.981 |
| eciRNA | F1 | 24 | 0.532 | 0.634 | +0.102 | 1.145 | 0.025 | 0.951 | 1×10⁻⁷ | 2.328 |
| EIciRNA | Precision | 24 | 0.629 | 0.740 | +0.112 | 1.107 | 0.044 | 0.802 | 4×10⁻⁴ | 0.870 |
| EIciRNA | Recall | 24 | 0.103 | 0.157 | +0.054 | 1.803 | −0.029 | 0.796 | 0.003 | 0.823 |
| EIciRNA | F1 | 24 | 0.168 | 0.243 | +0.075 | 1.659 | −0.036 | 0.784 | 8×10⁻⁴ | 0.857 |
| ciRNA | Precision | 24 | 0.754 | 0.885 | +0.131 | 0.595 | 0.437 | 0.953 | 6×10⁻⁵ | 1.000 |
| ciRNA | Recall | 24 | 0.114 | 0.352 | +0.238 | 2.854 | 0.026 | 0.925 | 2×10⁻⁵ | 1.161 |
| ciRNA | F1 | 24 | 0.182 | 0.442 | +0.259 | 2.028 | 0.072 | 0.913 | 2×10⁻⁵ | 1.430 |
| antisense | Precision | 21 | 0.663 | 0.837 | +0.173 | 0.475 | 0.522 | 0.406 | 1×10⁻⁶ | 0.521 |
| antisense | Recall | 24 | 0.135 | 0.356 | +0.221 | 2.035 | 0.081 | 0.629 | 1×10⁻⁷ | 1.015 |
| antisense | F1 | 24 | 0.216 | 0.450 | +0.234 | 1.344 | 0.160 | 0.578 | 1×10⁻⁷ | 1.083 |
| intergenic | Precision | 18 | 0.699 | 0.892 | +0.193 | 0.455 | 0.574 | 0.428 | 8×10⁻⁶ | 1.912 |
| intergenic | Recall | 24 | 0.089 | 0.217 | +0.128 | 2.119 | 0.029 | 0.920 | 1×10⁻⁴ | 1.169 |
| intergenic | F1 | 24 | 0.139 | 0.308 | +0.168 | 1.753 | 0.063 | 0.890 | 1×10⁻⁴ | 1.344 |

Every type's recall and F1 transfer well between the two benchmarks (R² 0.58–0.99): a caller that's relatively better or worse at a type on synthetic circRNAs stays that way on database circRNAs. For eciRNA, the database boost is close to the same fixed size for every caller. For EIciRNA, ciRNA, antisense, and intergenic, it isn't fixed: callers that were already doing well on synthetic circRNAs gain even more on database circRNAs, so the recall/F1 gap between strong and weak callers gets wider, not just shifted up equally for everyone. Precision behaves the opposite way for most types (R² 0.41–0.95): for eciRNA, ciRNA, antisense, and intergenic, weaker callers catch up instead of falling further behind, so precision differences between callers get smaller on database circRNAs. EIciRNA is the exception, it is the only type where precision differences widen. Every type shows a statistically significant improvement on database circRNAs (Wilcoxon p ≤ 3×10⁻³, Cohen's d 0.52–2.33), meaning they are easier to detect.

</details>

<details>
<summary><strong>Does caller ranking transfer?</strong></summary>

#### Does caller ranking transfer?

We averaged each caller's synthetic and database values across type and run first, ranked the 8 callers by that average on each side, then compared the two rankings with Kendall's tau-b:

| Metric | n callers | Kendall's tau-b | Spearman rho |
| --- | --- | --- | --- |
| Precision | 8 | 0.714 (p = 0.014) | 0.833 (p = 0.010) |
| Recall | 8 | 0.857 (p = 0.002) | 0.929 (p = 9×10⁻⁴) |
| F1 | 8 | 0.857 (p = 0.002) | 0.929 (p = 9×10⁻⁴) |

Ranking transfers significantly on every metric, including precision. By F1:

| Caller | F1 synthetic | F1 database | Rank synthetic | Rank database |
| --- | --- | --- | --- | --- |
| balanced_recall | 0.406 | 0.727 | 1 | 1 |
| discovery | 0.396 | 0.707 | 2 | 2 |
| CIRI-long | 0.319 | 0.470 | 3 | 4 |
| balanced_precision | 0.263 | 0.372 | 4 | 5 |
| circFL-seq | 0.251 | 0.538 | 5 | 3 |
| isoCirc | 0.188 | 0.339 | 6 | 6 |
| high_confidence | 0.102 | 0.106 | 7 | 7 |
| circNICK-lrs | 0.056 | 0.064 | 8 | 8 |

Six of the eight callers had the exact same rank on both benchmarks. circFL-seq and CIRI-long swapped places: circFL-seq ranked 5th on in silico circRNAs but 3rd on database ones, CIRI-long the reverse. This suggests that circFL-seq is stronger on known annotated circRNAs and worse on de novo. 

<img src="images/benchmark/slopegraph_database_human.png" width="100%"/>

</details>

<details>
<summary><strong>Per tool/mode, by metric</strong></summary>

#### Per tool/mode, by metric (n=15 paired points per row: the 5 circRNA types × 3 runs, except circNICK-lrs precision n=12 and high_confidence precision n=9)

| Tool | Metric | n | Synthetic | Database | Difference | Slope | Intercept | R² | Wilcoxon p | Cohen's d |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 15 | 0.687 | 0.865 | +0.178 | 0.286 | 0.669 | 0.217 | 6×10⁻⁵ | 1.664 |
| CIRI-long | Recall | 15 | 0.240 | 0.341 | +0.100 | 0.823 | 0.143 | 0.782 | 6×10⁻⁵ | 1.166 |
| CIRI-long | F1 | 15 | 0.319 | 0.470 | +0.151 | 1.045 | 0.136 | 0.750 | 6×10⁻⁵ | 1.707 |
| isoCirc | Precision | 15 | 0.877 | 0.945 | +0.068 | 0.460 | 0.541 | 0.783 | 6×10⁻⁵ | 1.405 |
| isoCirc | Recall | 15 | 0.121 | 0.217 | +0.096 | 0.769 | 0.124 | 0.874 | 6×10⁻⁵ | 1.882 |
| isoCirc | F1 | 15 | 0.188 | 0.339 | +0.152 | 0.644 | 0.218 | 0.822 | 6×10⁻⁵ | 1.699 |
| circNICK-lrs | Precision | 12 | 0.083 | 0.213 | +0.130 | −0.246 | 0.234 | 0.017 | 0.064 | 0.654 |
| circNICK-lrs | Recall | 15 | 0.056 | 0.056 | −0.000 | 1.115 | −0.007 | 0.980 | 0.530 | −0.008 |
| circNICK-lrs | F1 | 15 | 0.056 | 0.064 | +0.008 | 1.351 | −0.011 | 0.917 | 0.239 | 0.185 |
| circFL-seq | Precision | 15 | 0.624 | 0.922 | +0.298 | −0.075 | 0.969 | 0.146 | 6×10⁻⁵ | 0.818 |
| circFL-seq | Recall | 15 | 0.174 | 0.392 | +0.218 | 0.527 | 0.300 | 0.434 | 6×10⁻⁵ | 1.614 |
| circFL-seq | F1 | 15 | 0.251 | 0.538 | +0.287 | 0.331 | 0.455 | 0.229 | 6×10⁻⁵ | 1.513 |
| high_confidence | Precision | 9 | 0.900 | 0.939 | +0.039 | 0.655 | 0.349 | 0.931 | 0.024 | 1.302 |
| high_confidence | Recall | 15 | 0.065 | 0.068 | +0.003 | 1.089 | −0.003 | 0.996 | 0.762 | 0.241 |
| high_confidence | F1 | 15 | 0.102 | 0.106 | +0.004 | 1.085 | −0.005 | 0.993 | 0.762 | 0.188 |
| balanced_precision | Precision | 15 | 0.888 | 0.955 | +0.067 | 0.349 | 0.645 | 0.668 | 6×10⁻⁵ | 1.424 |
| balanced_precision | Recall | 15 | 0.186 | 0.254 | +0.069 | 0.911 | 0.085 | 0.951 | 6×10⁻⁵ | 1.557 |
| balanced_precision | F1 | 15 | 0.263 | 0.372 | +0.109 | 0.852 | 0.148 | 0.898 | 6×10⁻⁵ | 1.577 |
| balanced_recall | Precision | 15 | 0.690 | 0.895 | +0.205 | 0.337 | 0.663 | 0.950 | 6×10⁻⁵ | 2.353 |
| balanced_recall | Recall | 15 | 0.323 | 0.636 | +0.313 | 0.568 | 0.452 | 0.317 | 6×10⁻⁵ | 1.616 |
| balanced_recall | F1 | 15 | 0.406 | 0.727 | +0.321 | 0.772 | 0.414 | 0.452 | 6×10⁻⁵ | 2.501 |
| discovery | Precision | 15 | 0.650 | 0.813 | +0.163 | 0.716 | 0.348 | 0.850 | 6×10⁻⁵ | 2.336 |
| discovery | Recall | 15 | 0.331 | 0.647 | +0.316 | 0.581 | 0.455 | 0.336 | 6×10⁻⁵ | 1.629 |
| discovery | F1 | 15 | 0.396 | 0.707 | +0.311 | 0.873 | 0.361 | 0.389 | 6×10⁻⁵ | 2.274 |

isoCirc, high_confidence, and balanced_precision are consistent on all three metrics (R² ≥ 0.67). balanced_recall and discovery are strong on precision (R² ≥ 0.85) but less consistent on recall and F1 (R² 0.32–0.45, still significant). CIRI-long and circNICK-lrs are the opposite: precision doesn't track at all, but recall and F1 do (R² ≥ 0.75). circFL-seq is the least consistent of all eight: only its recall fit is significant, and even that is weak (R² = 0.43). By F1, high_confidence and circNICK-lrs perform almost identically on both benchmarks (difference ≤ 0.01), while balanced_recall, discovery, and circFL-seq score much higher on database circRNAs (+0.29 to +0.32); CIRI-long, isoCirc, and balanced_precision fall in between (+0.11 to +0.15).

| Tool | Outcome | Strata (types) | OR (MH) | 95% CI | MH p | Breslow-Day p |
| --- | --- | --- | --- | --- | --- | --- |
| CIRI-long | Precision | 5 | 3.90 | 3.62–4.20 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| CIRI-long | Recall | 5 | 1.73 | 1.67–1.81 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| isoCirc | Precision | 5 | 1.95 | 1.69–2.25 | <10⁻³⁰⁰ | 0.032 |
| isoCirc | Recall | 5 | 2.07 | 1.97–2.17 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circNICK-lrs | Precision | 5 | 1.54 | 1.43–1.67 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circNICK-lrs | Recall | 5 | 1.00 | 0.92–1.08 | 0.973 | <10⁻³⁰⁰ |
| circFL-seq | Precision | 5 | 2.91 | 2.65–3.19 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| circFL-seq | Recall | 5 | 3.34 | 3.20–3.49 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| high_confidence | Precision | 5 | 2.15 | 1.73–2.68 | 4×10⁻¹² | 0.247 |
| high_confidence | Recall | 5 | 1.07 | 0.99–1.15 | 0.082 | 1×10⁻¹⁶ |
| balanced_precision | Precision | 5 | 2.83 | 2.51–3.19 | <10⁻³⁰⁰ | 0.008 |
| balanced_precision | Recall | 5 | 1.60 | 1.52–1.67 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| balanced_recall | Precision | 5 | 3.98 | 3.74–4.24 | <10⁻³⁰⁰ | 7×10⁻⁵ |
| balanced_recall | Recall | 5 | 4.21 | 4.05–4.38 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Precision | 5 | 2.44 | 2.32–2.57 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |
| discovery | Recall | 5 | 4.33 | 4.17–4.51 | <10⁻³⁰⁰ | <10⁻³⁰⁰ |

Every caller has significantly higher odds of a correct call (precision) and of recovering a ground-truth circRNA (recall) on database circRNAs than on synthetic ones (OR 1.5–4.3), except circNICK-lrs's recall, which shows no real difference (OR ≈ 1.00, matching its flat recall line in the slopegraph). But that boost isn't the same size across circRNA types: of the sixteen caller/outcome combinations tested, only high_confidence's precision is uniform across types (Breslow-Day p = 0.25); every other combination's boost depends on which circRNA type it is.

</details>

**Conclusion.** Correlation is good: aggregate R² is 0.61–0.69, and most types and callers individually fit a significant line. Caller ranking transfers: Kendall's tau-b is 0.71–0.86 on every metric, including precision, and six of the eight callers land on the same F1 rank on both benchmarks, with only circFL-seq and CIRI-long swapping places with each other. But the size of the database boost isn't uniform: it depends on circRNA type for fifteen of the sixteen caller/outcome combinations tested (high_confidence's precision is the one exception). So the overall pattern and ranking are preserved well between de novo and database calling, but the exact performance improvement on annotated circRNAs is not uniform, and depends on both the circRNA type and the caller.

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

The annotated vs de novo distance (≈19) is not disproportionately large next to other, it's in the same order of magnitude. That supports the earlier hypothesis that de novo synthetic circRNAs are a valid comparator: they're offset from database circRNAs in absolute difficulty (expected as per reasons already documented above), by a margin similar in scale to other differences (species, protocol) while caller ranking and per-type structure transfer well between the two (see [Does the synthetic benchmark conclusions transfer to the database-derived dataset?](#does-the-synthetic-benchmark-conclusions-transfer-to-the-database-derived-dataset)).

---

For output file formats see [docs/output.md](output.md). For parameters see [docs/usage.md](usage.md). For merge/scoring algorithms see [docs/methods.md](methods.md).
