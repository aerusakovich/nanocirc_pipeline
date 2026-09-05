# Example output

Real output from a run of the pipeline's own bundled chr21 test dataset
(`assets/test_data/human_chr21/`), with cross-run merging and quantification
both enabled, so it also covers the parts the minimal 1-sample test in
[docs/quickstart.md](../quickstart.md#8-test-your-setup-first) does not.
Reproducible with the exact command in
[docs/quickstart.md, section 9](../quickstart.md#9-full-example-cross-run-merge-and-quantification).

Two samples (`chr21_smoke_a`, `chr21_smoke_b`), same input reads, treated as
2 independent runs of one group (`smoke`) so cross-run merging has something
to merge. All 4 detection tools and all 3 nanocirc confidence tiers
(Discovery, Balanced, High-Precision) ran.

## Files

- `crossrun/{discovery,balanced,high_confidence}/`
  - `smoke_<tier>_crossrun_clean.tsv`: the filtered, wet-lab-friendly output
    per tier, one row per merged locus.
  - `smoke_<tier>_crossrun.bed12`: same loci, BED12.
- `quantify/`
  - `chr21_smoke_<a|b>_quant_final_counts.tsv`: per-sample, per-locus
    remap-based read counts (`--run_quantify true`).
- `deseq2/`
  - `deseq2_counts_<tier>.tsv`, `deseq2_coldata_<tier>.tsv`,
    `deseq2_features_<tier>.tsv`: DESeq2-ready matrix, sample metadata, and
    locus metadata, per tier.
- `multiqc_report.html`: aggregated QC report.
- `execution_report.html`, `pipeline_dag.html`: Nextflow's own resource
  usage report and pipeline DAG for this exact run.

See [docs/output.md](../output.md) for what each column means.
