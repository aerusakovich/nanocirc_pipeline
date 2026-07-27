# nanocirc: Benchmark Results and Recommendations

This benchmark is adapted from our own earlier benchmark study, Rusakovich et al., "Benchmarking circRNA detection tools from long-read sequencing using a data-driven and flexible simulation framework" (Peer Community Journal, 2026; see [CITATIONS.md](../CITATIONS.md) for the full reference), modified and extended based on our newer findings. 

This page summarizes the benchmark behind nanocirc's four confidence tiers and gives practical guidance on which one to use. Full methodology with the new simulator and benchmark scripts will be published to the benchmark git https://gitlab.com/bingo-igdr/nano-circ as v3. Results below are from simulated human and mouse datasets, 3 runs each.

## Which tier should I use?

Each tier trades recall against precision differently. There is no single best tier, only the right one for your question.

| Tier | What it optimizes for | Trade-off |
| ---- | ---------------------- | --------- |
| **`discovery`** | Recall. Keeps almost every circRNA any tool found. | Highest false-positive rate: computational artefacts and incorrectly predicted structures will pass alongside real calls. |
| **`high_confidence`** | Precision. Only calls with strong cross-tool agreement on both BSJ and structure or isoCirc calls that met defined read threshold survive. | Sacrifices recall to get the structure right: false positives are minimal, but the reduced set will miss real biology so it's not fit for exploratory studies. |
| **`balanced_recall`** | A middle ground leaning recall. Cuts most of `discovery`'s false-positive rate while keeping reasonable recall. | Good default for exploratory analyses where getting the exact internal exon structure right, or an occasional false positive, isn't critical to the biological question. |
| **`balanced_precision`** | A middle ground leaning precision. Keeps recall relatively high while pushing precision further than `balanced_recall`. | Best for precision-leaning analyses that need to avoid false positives but still want decent recall, not just the aggressively small `high_confidence` set. |

> **Recommendation:** use one of the two **balanced** tiers for most applications. Reserve `discovery` for exploratory passes over a dataset, and `high_confidence` for when you specifically need a small, near-certain set.

## Quantification advice

Quantification is its own remap step, independent of which tool originally called a locus: the reproducibility and expression correlation plots below show read counts agree well across tools and tiers, and against ground truth, so switching tier for a re-analysis doesn't mean changing to a different counting method.

Please note that for differential expression (read counts, DESeq2), `discovery`'s false positives are not easy-to-spot near-zero-read artefacts that low-count filtering could remove (median ~16 reads vs ~16-23 for true positives in our benchmark), so they will be present in the count matrix as real-looking features. `high_confidence`, on the other end, drops enough real calls to hurt power on lower-expressed circRNAs. This is why we recommend the Balanced modes for expression-based analyses.

For a multi-sample design, turn on `--run_crossrun_merge` alongside `--run_quantify` so every sample in a `group` is quantified against one shared catalog. Otherwise each sample gets its own catalog and rows don't line up across samples in the [DESeq2 matrix](output.md#deseq2-count-matrix).

## Benchmark results

### Overall precision/recall

<img src="images/benchmark/pr_scatter_human.png" width="100%"/>
<img src="images/benchmark/pr_scatter_mouse.png" width="100%"/>

Precision/recall for every tool and every nanocirc tier, at three structural matching strictness levels (exon-based, relaxed BSJ, strict BSJ). The four tiers trace out the trade-off described above: `discovery` sits at high recall/lower precision, `high_confidence` at the opposite corner, and the two balanced tiers in between.

### Ground truth coverage

<img src="images/benchmark/gt_coverage_stacked_human.png" width="100%"/>
<img src="images/benchmark/gt_coverage_stacked_mouse.png" width="100%"/>

Green (recovered) plus grey (missed) always sum to exactly 100% of ground truth, so that split directly shows a caller's recall. Red stacks separately on top of the 100% line as false calls the caller made scaled the same way. 

### What are these false positives?

We traced every FP-supporting read back to the simulator's own ground truth (single-tool discovery exon-based FPs, all 3 human runs, mean % ± SD across runs):

| FP read is... | CIRI-long | CircFL-seq | IsoCirc | CircNick-LRS |
| --- | --- | --- | --- | --- |
| Real circ read, correct boundary, wrong strand (antisense called as sense) | 86.0 ± 0.6 | 79.4 ± 1.7 | 18.8 ± 7.6 | 9.9 ± 0.6 |
| Real circ read, correct strand, pure boundary noise | 0.7 ± 0.4 | 0.0 ± 0.0 | 0.0 ± 0.0 | 43.8 ± 0.7 |
| Real circ read, but circRNA reported at a different/distant locus | 13.4 ± 0.7 | 20.6 ± 1.7 | 81.2 ± 7.6 | 46.3 ± 0.3 |
| Not from a circ read at all (linear read miscalled) | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.1 ± 0.1 |
| FP locus overlaps a genomic repeat | 19.8 ± 0.0 | 39.9 ± 0.7** | **55.2 ± 8.4** | **77.9 ± 0.2** |
| n FP loci (pooled, 3 runs) | 1,534 | 559 | 21 | 2,754 |
| n supporting reads (pooled, 3 runs) | 32,593 | 3,687 | 25 | 11,025 |


Many FP come from `antisense`-type circRNA, where the simulator deliberately sets the circRNA's strand opposite to its host gene. Simulator also mimics cDNA based RCRT-sequencing, meaning it simulates both strand sequences: one as an actual circRNA read and the other as cDNA read. The tools correctly find the genomic position of circRNA but select the host gene's (wrong) strand instead of flipping it to the antisense one (correct). Since strand agreement is required by all three matching strategies (exon-based, relaxed BSJ, strict BSJ), no boundary-tolerance setting can rescue these FP and we didn't find an easy way to override this gene-strand bias without introducing new FP elsewhere. These realistic strand misses could also happen on real data and we suggest to keep that in mind as you move onto wet-lab validation. 

**Practical implication**: most `discovery` false positives still point at the right gene/circRNA locus, just with an imprecise boundary or wrong strand. That makes them usable for gene- or locus-level expression signal analysis, but not for isoform-level analysis or primer/probe design, both of which need to get the exact strand/boundary right.

### What discovery still misses?

We also characterized every ground-truth circRNA against whether any tool recovered it at all (mean % ± SD across 3 human runs; missed + recovered = 100% of the GT catalog, ~8,700 circRNAs/run):

| | Missed (FN) | Recovered |
| --- | --- | --- |
| Share of GT catalog | 63.1 ± 0.1 | 36.9 ± 0.1 |
| Median mature length (bp) | 1,057 ± 12 | 291 ± 1 |
| Mean rolling-circle copy number | 7.6 ± 0.1 | 8.8 ± 0.0 |
| Median simulated read count | 34 ± 0 | 99 ± 1 |
| Canonical GT/AG splice site at BSJ | 63.0 ± 0.0 | 81.0 ± 0.0 |
| No gene overlap at all (intergenic) | 20.0 ± 0.0 | 9.0 ± 0.0 |
| Falls in a gene with ≥5 isoforms | 73.0 ± 0.0 | 84.0 ± 0.0 |
| Overlaps a genomic repeat | 51.0 ± 1.0 | 58.0 ± 1.0 |

Recall per biotype:

| Biotype | Missed | Recovered |
| --- | --- | --- |
| `eciRNA` (pure exonic) | 20.7 ± 0.2 | 79.3 ± 0.2 |
| `antisense` | 68.9 ± 0.1 | 31.1 ± 0.1 |
| `ciRNA` (purely intronic) | 74.5 ± 0.2 | 25.5 ± 0.2 |
| `EIciRNA` (exon+intron) | 76.4 ± 0.1 | 23.6 ± 0.1 |
| `intergenic` | 81.3 ± 0.2 | 18.7 ± 0.2 |

Missed circRNAs are generally longer (~3.6x), lower-expressed (~3x fewer reads), and use non-canonical splice sites more often. Rolling-circle copy number is also slightly lower, meaning weaker RCRT/RCA tandem-repeat evidence for the tools that rely on it. `eciRNA` seems to be the easier type to recover; every other biotype, especially without a clean match to an annotated exon, is more likely to be missed. Repeat overlap isn't a meaningful driver of circRNA discovery. 

### circRNA type classification

<img src="images/benchmark/type_classification_human_exon.png" width="100%"/>
<img src="images/benchmark/type_classification_human_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_human_strict_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_exon.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_relaxed_bsj.png" width="100%"/>
<img src="images/benchmark/type_classification_mouse_strict_bsj.png" width="100%"/>

Same TP/FP-vs-GT breakdown as the GT coverage plot above, split by circRNA biotype (`eciRNA`, `EIciRNA`, `ciRNA`, `antisense`, `intergenic`, see [docs/output.md](output.md#circrna-type-classification)) instead of one bar per caller.

### Expression quantification and reproducibility

<img src="images/benchmark/expression_correlation_human.png" width="100%"/>
<img src="images/benchmark/expression_correlation_mouse.png" width="100%"/>

Read-count correlation against ground truth (right panel) and cross-method agreement (left panel), per tool/tier. Each cell shows three metrics:

- **Pearson r**: linear correlation between predicted and true counts. Sensitive to whether counts scale proportionally, not just in the same direction.
- **Spearman rho**: rank correlation. Only checks whether higher-expressed circRNAs get higher counts, regardless of exact scale. More resistant to outliers and non-linear (but still monotonic) relationships.
- **R²** (Pearson r squared): how much of the variance in true expression is captured by the predicted counts, on a 0-1 scale.

For differential expression, Spearman rho matters to rank circRNAs correctly relative to each other; Pearson r/R² matter to ensure the fold-changes DESeq2 estimates between conditions are quantitatively trustworthy, not just directionally right.

<img src="images/benchmark/sensitivity_by_expression_human.png" width="100%"/>
<img src="images/benchmark/sensitivity_by_expression_mouse.png" width="100%"/>

Detection sensitivity broken down by ground-truth expression tertile: low-expressed circRNAs are structurally harder to detect correctly than high-expressed ones, for every tool.

<img src="images/benchmark/reproducibility_pairs_human.png" width="100%"/>
<img src="images/benchmark/reproducibility_pairs_mouse.png" width="100%"/>

Pairwise read-count agreement across every tool/tier combination, ground truth included.

---

For output file formats see [docs/output.md](output.md). For parameters see [docs/usage.md](usage.md).
