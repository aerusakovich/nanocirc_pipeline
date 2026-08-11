# nanocirc: Methods

This page describes how nanocirc's merging, confidence scoring, cross-run merge, quantification, and type classification actually work. For command-line parameters, see [docs/usage.md](usage.md). For output file locations and column schemas, see [docs/output.md](output.md). For benchmark results and tier recommendations, see [docs/benchmark.md](benchmark.md).

## Confidence scoring

Each merged circRNA is scored on two **independent** axes:

- **`bsj_consensus`**: what fraction of active tools agreed on this back-splice junction?
- **`isoform_consensus`**: how many tools confirmed a matching exon structure at that BSJ?

A circRNA can have a well-supported BSJ (`High` `bsj_consensus`) but an uncertain isoform structure (`Low` `isoform_consensus`), or the other way round.

### Scoring bins

`bsj_score`/`isoform_score` are a direct function of the supporting-tool count, independent of how many tools were active (`n_active`):

| Supporting tool count | Score | Consensus |
| ---------------------- | ----- | --------- |
| 1                       | 1     | Low       |
| 2                       | 2     | Medium    |
| 3 or more               | 4     | High      |

A given score means the same amount of evidence in every run, regardless of `n_active`. This allows the confidence score to be used as a stable metric across runs with different tool counts.

### Scoring examples

Verified against `count_to_score()`/`score_to_cat()` in `bin/add_isoform_confidence.py`.

#### 4-tool run

| Scenario                               | BSJ tools  | Isoform tools | bsj score | bsj_consensus | iso score | isoform_consensus |
| --------------------------------------- | ---------- | -------------- | --------- | -------------- | --------- | ------------------ |
| All 4 agree, full isoform match        | 4          | 4              | 4         | **High**       | 4         | **High**            |
| 3 tools agree, good isoform            | 3          | 3              | 4         | **High**       | 4         | **High**            |
| All 4 BSJ, but no isoform confirmation | 4          | 0              | 4         | **High**       | 1         | **Low**             |
| 2 tools agree, no isoform confirmation | 2          | 0              | 2         | **Medium**     | 1         | **Low**             |
| Only 1 tool detects this circRNA       | 1          | 1              | 1         | **Low**        | 1         | **Low**             |

#### 3-tool run

| Scenario                               | BSJ tools | Isoform tools | bsj score | bsj_consensus | iso score | isoform_consensus |
| --------------------------------------- | --------- | -------------- | --------- | -------------- | --------- | ------------------ |
| All 3 agree, full isoform match        | 3         | 3              | 4         | **High**       | 4         | **High**            |
| 2 tools agree, good isoform agreement           | 2         | 2              | 2         | **Medium**     | 2         | **Medium**          |
| All 3 BSJ, but no isoform agreement | 3         | 0              | 4         | **High**       | 1         | **Low**             |
| Only 1 tool                            | 1         | 0              | 1         | **Low**        | 1         | **Low**             |

#### 2-tool run

| Scenario                        | BSJ tools | Isoform tools | bsj score | bsj_consensus | iso score | isoform_consensus |
| --------------------------------- | --------- | -------------- | --------- | -------------- | --------- | ------------------ |
| Both agree, full isoform match  | 2         | 2              | 2         | **Medium**     | 2         | **Medium**          |
| Both agree, no isoform          | 2         | 0              | 2         | **Medium**     | 1         | **Low**             |
| Only 1 tool detects             | 1         | 0              | 1         | **Low**        | 1         | **Low**             |

> [!WARNING]
> `High` requires 3+ agreeing tools regardless of how many were active, so `High` is never reachable in a 1- or 2-tool run, even if every active tool agrees. Always consider the raw `bsj_confidence`/`isoform_confidence` tool counts alongside the consensus labels, not just the label itself.

## Merge algorithms

All tools within a relaxed-BSJ group (coordinates within `--circrna_bsj_tolerance` bp) are treated as candidates for the same circRNA. A merge algorithm selects the representative BSJ and exon structure from those candidates.

**BSJ selection** (all modes except `priority`): majority vote across all tools; ties broken by tool priority IsoCirc > CircFL-seq > CIRI-long > CircNick-LRS.

### Default merge mode: `consensus_hybrid`

The pipeline uses `consensus_hybrid` for all four confidence tiers: `discovery`, `balanced_precision`, `balanced_recall`, and `high_confidence`.

| Property | `consensus_hybrid` |
| -------- | ------------------- |
| BSJ | Majority vote across all tools |
| Structure vote participants | Tools sharing the **exact winning BSJ** only |
| Coordinate comparison | Absolute genomic coords (boundaries within `--circrna_bsj_tolerance` bp) |
| Rebasing of minority-BSJ tools | No: tools with a different BSJ do not contribute to the structure vote |
| Minority-BSJ tool handling | Emitted as separate isoform entries at their own BSJ coordinates |
| Structure tie-break priority | IsoCirc > CircFL-seq > CIRI-long > CircNick-LRS |

This selects the most-supported exon structure among tools that agree on the BSJ, without shifting coordinates from tools that chose a slightly different junction. Minority-BSJ isoforms are preserved in the output rather than discarded.

**Multi-isoform recovery**: `isocirc`, `circfl-seq`, and `ciri-long` each often correctly resolve more than one real isoform at a single BSJ (different reads supporting different structures). Rather than collapsing each tool to its single best call before voting, `consensus_hybrid` runs the vote above to pick the winning BSJ and its primary structure, then adds back every other isoform call from those three tools at the winning BSJ that doesn't already match an existing entry. `circnick-lrs` collapses to its single best call, since its own isoform-level calls are unreliable (near-zero true-isoform recovery even under a loose reciprocal-overlap threshold). 

### Additional modes (`--run_benchmark_modes` only)

Three further algorithms were tested and outperformed by `consensus_hybrid` in benchmarking, and are neither computed nor published by default (skipped entirely unless `--run_benchmark_modes true`, and never computed in cross-run mode). Kept for internal benchmarking history and as diagnostic baselines, not maintained for routine use.

| Mode | BSJ selection | Structure vote participants | Minority-BSJ rebasing | Minority-BSJ handling |
| ---- | -------------- | ---------------------------- | ----------------------- | ----------------------- |
| `consensus` | Majority vote | Tools sharing the exact winning BSJ only | No | Separate isoforms at own coords |
| `consensus_xstruct` | Majority vote | All tools in group | Yes, shifted to winning BSJ | Folded into structure vote |
| `priority` | Highest-priority tool | Highest-priority tool | No | Separate isoforms at own coords |

- `consensus` is structurally the same vote as `consensus_hybrid`, but compares structure by string equality rather than coordinate similarity: two tools reporting the same exon structure with a 1bp boundary difference are treated as different isoforms. It also does not get the multi-isoform recovery pass described above.
- `consensus_xstruct` includes minority-BSJ tools in the structure vote by rebasing their exon coordinates to the winning BSJ. This can incorporate more structural information but may introduce coordinate imprecision when the BSJ offset is non-trivial.
- `priority` skips voting entirely: BSJ and structure come unconditionally from the single highest-priority tool present.

### Confidence tiers

After merging, each entry is scored on the two independent axes above, then one of four filters is applied:

| Tier | Merge algorithm | Filter | Rule | Axes retained |
| ---- | ----------------- | ------ | ---- | --------------- |
| `discovery` | `consensus_hybrid` | none | Keep all entries | any |
| `balanced_precision` | `consensus_hybrid` | `isocirc_only` | Drop entries where either axis is Low, unless the source is IsoCirc | ≥ Medium on both, or Low from IsoCirc |
| `balanced_recall` | `consensus_hybrid` | `trusted_only` | Drop entries where either axis is Low, unless the source is a trusted tool (`--circrna_trusted_tools`, default CIRI-long/IsoCirc/CircFL-seq) | ≥ Medium on both, or Low from a trusted tool |
| `high_confidence` | `consensus_hybrid` | `high_only_isocirc` | Keep only entries where both axes are High, unless Low and the source is IsoCirc | High on both, or Low from IsoCirc |

If `--run_quantify true` is set, `discovery`/`balanced_recall`/`high_confidence` get a further post-quantification confidence filter on top of this table, see [Quantification](#quantification) below. It does not run on `balanced_precision`: its Low-from-IsoCirc entries pass as-is, favoring recall, whereas `high_confidence` additionally requires read support there since that tier prioritizes precision over recall.

### Why we chose vote-based merging and not simple union/intersection approaches? 

A natural alternative to nanocirc's tiers is to just deduplicate each tool's raw calls (relaxed-BSJ grouping) and take the union (any 1 of 4 tools) or an N-of-4 intersection. Scored the same way (exon-based matching, `final_benchmark` human/mouse, pooled TP/FP/FN across 3 runs each) against the matched nanocirc tier for each union/intersection rule:

**Sensitivity: union vs `discovery`**

| | Human | Mouse |
| --- | --- | --- |
| Naive union (any 1 of 4 tools, deduplicated) | P 0.546 / R 0.381 / F1 0.449 | P 0.602 / R 0.527 / F1 0.562 |
| `discovery` | P 0.577 / R 0.392 / F1 0.467 | P 0.665 / R 0.556 / F1 0.606 |

**Precision: intersection vs `balanced_precision`/`high_confidence`**

| | Human | Mouse |
| --- | --- | --- |
| Naive 3-of-4 intersection | P 0.892 / R 0.107 / F1 0.192 | P 0.980 / R 0.159 / F1 0.274 |
| `balanced_precision` | P 0.850 / R 0.202 / F1 0.327 | P 0.942 / R 0.361 / F1 0.522 |
| Naive 4-of-4 intersection | P 0.911 / R 0.050 / F1 0.095 | P 0.997 / R 0.086 / F1 0.158 |
| `high_confidence` | P 0.910 / R 0.095 / F1 0.171 | P 0.986 / R 0.196 / F1 0.328 |

**F1: 2-of-4 intersection vs `balanced_recall`**

| | Human | Mouse |
| --- | --- | --- |
| Naive 2-of-4 intersection | P 0.827 / R 0.186 / F1 0.304 | P 0.910 / R 0.285 / F1 0.434 |
| `balanced_recall` | P 0.740 / R 0.360 / F1 0.484 | P 0.849 / R 0.522 / F1 0.646 |

At every matched threshold, nanocirc's tier has higher recall and higher F1 than the basic operation count, in both species. nanocirc's tiers separate the BSJ-agreement decision from the structure-agreement decision (a locus can qualify without every supporting tool agreeing on the exact exon structure) and recover a tool's own additional correctly-resolved isoforms at the winning BSJ (see [multi-isoform recovery](#default-merge-mode-consensus_hybrid) above). That is the main reason nanocirc's consensus has a better precision/recall balance than simple set operations on the same 4 tools' raw calls.

## Cross-run merge

When `--run_crossrun_merge true` is set and the samplesheet has a `group` column, all samples sharing a group are merged after per-sample analysis. Cross-run merge runs once per tier: each tier's crossrun input is that tier's own already-filtered per-sample output (e.g. `balanced_recall`'s crossrun merge starts from every run's own `trusted_only`-filtered `balanced_recall` catalog, not from `discovery`). Each run is then treated as an independent caller within that tier, using the same `consensus_hybrid` vote and confidence scoring used across tools within one sample, treating samples as tools.

Each candidate locus is scored with two separate votes, run independently:

1. **BSJ vote**: counts how many different runs support the exact back-splice junction position, ignoring exon structure. Decides which tier a locus qualifies for.
2. **Structure vote**: restricted to records at the winning BSJ position only, groups them by exon structure and ranks groups by total tool agreement. Decides the final isoform structure reported for the locus.

Splitting these two votes stops a locus's correct BSJ call from being discarded just because a different run reported a slightly different exon structure at the same junction.

A structure seen in only 1 run is dropped from the structure vote unless that run's own tool agreement meets `--crossrun_min_tool_agreement` (default `2`); circRNAs with evidence in 2 or more runs are never dropped by this check. This removes weak, single-run calls without changing which tier the underlying BSJ qualifies for, and cut cross-run false positives by about 34% in testing, for a small recall cost. Minority BSJ positions (backed by fewer runs than the winning position) each keep their own separate isoform entry rather than being merged.

| Tier | What it takes to be retained |
| ----- | ------------------------------ |
| `discovery` | Detected by at least 1 tool in at least 1 run. Maximum sensitivity. |
| `balanced_precision` | Multiple tools agreed within a run (or a Low-confidence IsoCirc-only call), and multiple runs support it. |
| `balanced_recall` | Same run-count threshold as `balanced_precision`; built from each run's own `trusted_only`-filtered `balanced_recall` catalog, so it allows per-sample tier's Low-confidence calls from trusted tools that `balanced_precision`'s per-sample `isocirc_only` filter would drop. |
| `high_confidence` | Built from each run's own `high_only_isocirc`-filtered `high_confidence` catalog (most/all tools agreed within each run), and most/all runs support it. Maximum precision. |

**Count thresholds**, where `n` is the number of runs in the group, verified against `circrna_crossrun_merge.nf`:

| Tier | Minimum runs required |
| ----- | ------------------------ |
| `discovery` | ≥ 1 (all circRNAs retained) |
| `balanced_precision` | ≥ max(2, ceil(0.25 × n)) |
| `balanced_recall` | ≥ max(2, ceil(0.25 × n)) |
| `high_confidence` | ≥ ceil(0.75 × n) |

## Quantification

When `--run_quantify true` is set, each sample's reads are remapped against synthetic tandem-duplicated circle references built from a `discovery` (unfiltered, hybrid-consensus) locus catalog:

1. **Tier1**: chunked remap-and-classify pass against reference chunks (`--quant_chunk_size` loci per chunk), plus an overlap-cluster rescue pass for near-duplicate loci (`--quant_locus_dedup_tolerance` bp).
2. **Tier2**: targeted low-coverage rescue for loci tier1 handled poorly (flagged via `--quant_min_old_tool_count`).
3. **Tier3**: gene-family/repeat-cluster rescue for loci that remain ambiguous after tier2.

**Crossrun approach** toggled by `--run_crossrun_merge`. If `off`, each sample is quantified against its own per-sample `discovery` catalog. If `on`, every sample sharing a `group` is quantified against one shared catalog built from that group's crossrun `discovery` output, so circRNAs that weren't called by a tool in this run can still be recovered and quantified if they were called in another run of this group. The circle references, chunk references, and gene-family similarity database are built once per group, not once per sample.

### Post-quantification confidence filter

After quantification, `discovery`/`balanced_recall`/`high_confidence` get an additional filter based on the quantified read count (`--circrna_confident_min_reads`, default `2`):

- `discovery`/`balanced_recall`: drops loci only CircNick-LRS called (`supporting_tools` is just `circnick`) with a read count at or below the threshold.
- `high_confidence`: drops IsoCirc-only calls (its merge-time `high_only_isocirc` exception) at or below the same threshold, since this tier favors precision over recall.
- `balanced_precision` is not filtered this way: its IsoCirc-only entries (same merge-time exception) pass through unconditionally, since this tier favors recall.

## circRNA type classification

Types are assigned by intersecting BSJ coordinates against gene and exon BED files derived from the GTF:

| Type | Definition |
| ----- | ----------- |
| `eciRNA` | Same-strand gene overlap; fully covered by exonic regions (purely exonic) |
| `EIciRNA` | Same-strand gene overlap; overlaps exons but retains intronic content (exon-intron circRNA) |
| `ciRNA` | Same-strand gene overlap; no exon overlap (purely intronic) |
| `antisense` | Overlaps a gene on the opposite strand only |
| `intergenic` | No overlap with any annotated gene on either strand |