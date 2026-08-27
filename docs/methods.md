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

The pipeline uses `consensus_hybrid` for all three confidence tiers: `discovery`, `balanced`, and `high_confidence`.

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

### Splice-motif strand recheck

When a reference FASTA and its `.fai` index are available (`--fasta`), every merged entry's called strand is checked against the genomic splice-site motif at each of its exon-exon junctions, including the outer BSJ: `+` strand expects donor/acceptor `GT`/`AG`, `-` strand expects the reverse-complement `CT`/`AC`. This compares only against the tool's own called strand, never against ground truth (unavailable on real data).

- Single-junction (single-exon) calls are never hard-flipped to the opposite strand: on real antisense-vs-sense-gene data, a lone 2bp motif check is confidently wrong about 10% of the time, too risky to act on alone. They can only be confirmed or downgraded to ambiguous (`.`).
- Calls with 2 or more junctions are flipped to the opposite strand only when every junction unanimously supports it and none support the called strand; any other mixed signal downgrades to ambiguous (`.`) instead of flipping.
- A call with no strand at all (`.`) is left alone: there is nothing to check.
- Skipped entirely, with every entry passed through unchecked, when the FASTA has no `.fai` next to it, since the recheck isn't required for the pipeline to run.

**Collision handling**: two entries can independently resolve to the same `(chrom, start, end, strand)` -- most commonly two different candidate groups both downgrading to the same ambiguous `.` coordinates -- which would otherwise produce duplicate `bsj_id`s, a uniqueness every downstream quantification script assumes holds. An entry that collides is pushed to ambiguous (`.`) first; if it still collides there, a `dup{N}` suffix (`dup2`, `dup3`, ...) is appended to its `bsj_id` to guarantee uniqueness.

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

After merging, each entry is scored on the two independent axes above, then one of three filters is applied:

| Tier | Merge algorithm | Filter | Rule | Axes retained |
| ---- | ----------------- | ------ | ---- | --------------- |
| `discovery` | `consensus_hybrid` | none | Keep all entries | any |
| `balanced` | `consensus_hybrid` | `trusted_only` | Drop entries where either axis is Low, unless the source is a trusted tool (`--circrna_trusted_tools`, default CIRI-long/IsoCirc/CircFL-seq) | ≥ Medium on both, or Low from a trusted tool |
| `high_confidence` | `consensus_hybrid` | `isocirc_only` | Drop entries where either axis is Low, unless the source is IsoCirc | ≥ Medium on both, or Low from IsoCirc |

If `--run_quantify true` is set, all three tiers get a further post-quantification confidence filter on top of this table, see [Quantification](#quantification) below. On `discovery`/`balanced` this drops loci where CircNick-LRS is the only supporting tool and read support is weak; on `high_confidence` it instead drops loci where IsoCirc is the only supporting tool and read support is weak, guarding the "Low from IsoCirc" exception opened by `isocirc_only` above.

### Why we chose vote-based merging and not simple union/intersection approaches? 

A natural alternative to nanocirc's tiers is to just deduplicate each tool's raw calls (relaxed-BSJ grouping) and take the union (any 1 of 4 tools) or an N-of-4 intersection. Scored the same way (exon-based matching, `final_benchmark` human/mouse, pooled TP/FP/FN across 3 runs each) against the matched nanocirc tier for each union/intersection rule:

**Sensitivity: union vs `discovery`**

| | Human | Mouse |
| --- | --- | --- |
| Naive union (any 1 of 4 tools, deduplicated) | P 0.496 / R 0.339 / F1 0.402 | P 0.537 / R 0.467 / F1 0.500 |
| `discovery` | P 0.676 / R 0.369 / F1 0.477 | P 0.671 / R 0.493 / F1 0.569 |

**Precision: intersection vs `high_confidence`**

| | Human | Mouse |
| --- | --- | --- |
| Naive 3-of-4 intersection | P 0.729 / R 0.089 / F1 0.159 | P 0.814 / R 0.137 / F1 0.234 |
| `high_confidence` | P 0.856 / R 0.196 / F1 0.319 | P 0.900 / R 0.287 / F1 0.435 |

**F1: 2-of-4 intersection vs `balanced`**

| | Human | Mouse |
| --- | --- | --- |
| Naive 2-of-4 intersection | P 0.720 / R 0.178 / F1 0.286 | P 0.798 / R 0.265 / F1 0.398 |
| `balanced` | P 0.765 / R 0.358 / F1 0.488 | P 0.790 / R 0.483 / F1 0.600 |

At every matched threshold, nanocirc's tier has higher recall and higher F1 than the basic set operation, in both species, and higher precision too except on mouse's 2-of-4-vs-`balanced` pair, where the naive intersection is marginally ahead (0.798 vs 0.790). One caller row can only claim one ground-truth entry: several tools independently calling near-identical structures at the same locus do not inflate a naive union or intersection's true positives, they compete for the same GT entry instead, and the losers count as false positives. This is why naive precision here is lower than a simple tool-agreement count would suggest. nanocirc's tiers separate the BSJ-agreement decision from the structure-agreement decision (a locus can qualify without every supporting tool agreeing on the exact exon structure) and recover a tool's own additional correctly-resolved isoforms at the winning BSJ (see [multi-isoform recovery](#default-merge-mode-consensus_hybrid) above). That is the main reason nanocirc's consensus has a better precision/recall balance than simple set operations on the same 4 tools' raw calls.

## Cross-run merge

When `--run_crossrun_merge true` is set and the samplesheet has a `group` column, all samples sharing a group are merged after per-sample analysis. Cross-run merge runs once per tier: each tier's crossrun input is that tier's own already-filtered per-sample output (e.g. `balanced`'s crossrun merge starts from every run's own `trusted_only`-filtered `balanced` catalog, not from `discovery`). Each run is then treated as an independent caller within that tier, using the same `consensus_hybrid` vote and confidence scoring used across tools within one sample, treating samples as tools.

Each candidate locus is scored with two separate votes, run independently:

1. **BSJ vote**: counts how many different runs support the exact back-splice junction position, ignoring exon structure. Decides which tier a locus qualifies for.
2. **Structure vote**: restricted to records at the winning BSJ position only, groups them by exon structure and ranks groups by total tool agreement. Decides the final isoform structure reported for the locus.

Splitting these two votes stops a locus's correct BSJ call from being discarded just because a different run reported a slightly different exon structure at the same junction.

A structure seen in only 1 run is dropped from the structure vote unless that run's own tool agreement meets `--crossrun_min_tool_agreement` (default `2`); circRNAs with evidence in 2 or more runs are never dropped by this check. This removes weak, single-run calls without changing which tier the underlying BSJ qualifies for, and cut cross-run false positives by about 34% in testing, for a small recall cost. Minority BSJ positions (backed by fewer runs than the winning position) each keep their own separate isoform entry rather than being merged.

| Tier | What it takes to be retained |
| ----- | ------------------------------ |
| `discovery` | Detected by at least 1 tool in at least 1 run. Maximum sensitivity. |
| `balanced` | Built from each run's own `trusted_only`-filtered `balanced` catalog (multiple tools agreed within a run, or a Low-confidence call from a trusted tool), and multiple runs support it. |
| `high_confidence` | Built from each run's own `isocirc_only`-filtered `high_confidence` catalog (multiple tools agreed within a run, or a Low-confidence IsoCirc-only call), and multiple runs support it. |

**Count thresholds**, where `n` is the number of runs in the group, verified against `circrna_crossrun_merge.nf`:

| Tier | Minimum runs required |
| ----- | ------------------------ |
| `discovery` | ≥ 1 (all circRNAs retained) |
| `balanced` | ≥ max(2, ceil(0.25 × n)) |
| `high_confidence` | ≥ max(2, ceil(0.25 × n)) |

## Quantification

When `--run_quantify true` is set, each sample's reads are remapped against synthetic tandem-duplicated circle references built from a `discovery` (unfiltered, hybrid-consensus) locus catalog:

1. **Tier1**: chunked remap-and-classify pass against reference chunks (`--quant_chunk_size` loci per chunk), plus an overlap-cluster rescue pass for near-duplicate loci (`--quant_locus_dedup_tolerance` bp).
2. **Tier2**: targeted low-coverage rescue for loci tier1 handled poorly (flagged via `--quant_min_old_tool_count`).
3. **Tier3**: gene-family/repeat-cluster rescue for loci that remain ambiguous after tier2.

**Tier2/tier3 read attribution is boundary-aware.** Within a rescue group (sibling isoforms sharing a BSJ, or gene-family/repeat paralogs), a read that qualifies against more than one candidate reference is immediately assigned to the reference with the greatest total aligned span (matched CIGAR blocks only, deletions and introns excluded so a large deletion isn't counted as coverage) only if it actually includes the sequence unique to that reference: a k-mer (21bp by default) absent from every other candidate in the group. A read whose aligned span is entirely inside sequence shared by two or more candidates carries no information to decide between them, so instead of winner-take-all approach by total span it is split with equal fractional weight across all qualifying candidates (RSEM/Salmon-style multi-mapper handling). This is important for nested isoform pairs, where one isoform is a structural subset of the other: most of the shorter isoform's reads map almost entirely inside the region the two share, so total-span scoring would assign them to whichever candidate happens to align longer.

**Crossrun approach** toggled by `--run_crossrun_merge`. If `off`, each sample is quantified against its own per-sample `discovery` catalog. If `on`, every sample sharing a `group` is quantified against one shared catalog built from that group's crossrun `discovery` output, so circRNAs that weren't called by a tool in this run can still be recovered and quantified if they were called in another run of this group. The circle references, chunk references, and gene-family similarity database are built once per group, not once per sample.

### Post-quantification confidence filter

After quantification, `discovery`/`balanced`/`high_confidence` each get an additional filter based on the quantified read count (`--circrna_confident_min_reads`, default `2`):

- `discovery`/`balanced`: drops loci where CircNick-LRS is the only supporting tool (`supporting_tools` is just `circnick`) and the read count is at or below the threshold.
- `high_confidence`: drops loci where IsoCirc is the only supporting tool (`supporting_tools` is just `isocirc`) and the read count is at or below the threshold, guarding the merge-time `isocirc_only` exception that let that Low-confidence, single-tool call into this tier in the first place.

## circRNA type classification

Types are assigned by intersecting BSJ coordinates against gene and exon BED files derived from the GTF:

| Type | Definition |
| ----- | ----------- |
| `eciRNA` | Same-strand gene overlap; fully covered by exonic regions (purely exonic) |
| `EIciRNA` | Same-strand gene overlap; overlaps exons but retains intronic content (exon-intron circRNA) |
| `ciRNA` | Same-strand gene overlap; no exon overlap (purely intronic) |
| `antisense` | Overlaps a gene on the opposite strand only |
| `intergenic` | No overlap with any annotated gene on either strand |