# nanocirc: Usage

## Introduction

**nanocirc** is a pipeline for detection and characterisation of circular RNAs (circRNAs) from long-read nanopore sequencing data. It runs up to four detection tools in parallel, converts their outputs to a unified BED12 format, and merges results with confidence scoring when two or more tools are active.

## Samplesheet input

You must provide a samplesheet CSV file with two columns: `sample` and `fastq`.

```csv title="samplesheet.csv"
sample,fastq
SAMPLE1,/path/to/sample1.fastq.gz
SAMPLE2,/path/to/sample2.fastq.gz
```

| Column   | Description                                                                          |
| -------- | ------------------------------------------------------------------------------------ |
| `sample` | Unique sample name. Cannot contain spaces.                                           |
| `fastq`  | Full path to a gzipped FASTQ file (`.fastq.gz` or `.fq.gz`). Must be nanopore data. |

Pass the samplesheet to the pipeline with:

```bash
--input '[path to samplesheet.csv]'
```

## Running the pipeline

### Minimal example

At minimum you need the samplesheet, a reference genome FASTA, a GTF annotation, and at least one detection tool enabled:

```bash
nextflow run aerusakovich/nanocirc_pipeline \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results/ \
    --fasta /path/to/genome.fa \
    --gtf   /path/to/annotation.gtf
```

### Full example (all four tools)

```bash
nextflow run aerusakovich/nanocirc_pipeline \
    -profile singularity \
    --input        samplesheet.csv \
    --outdir       results/ \
    --fasta        /path/to/genome.fa \
    --gtf          /path/to/annotation.gtf \
    --circrna_db   /path/to/circ_db.bed \
    --run_isocirc  true \
    --run_circfl   true \
    --run_cirilong true \
    --run_circnick true \
    --circnick_species mouse \
    -resume
```

> [!NOTE]
> `--circrna_db` is required when `--run_isocirc` or `--run_cirilong` is set.

### Resuming a run

Add `-resume` to any command to reuse cached results from a previous run:

```bash
nextflow run aerusakovich/nanocirc_pipeline ... -resume
```

### Params file

For repeated runs with the same settings, use a params YAML file:

```bash
nextflow run aerusakovich/nanocirc_pipeline -profile singularity -params-file params.yaml
```

```yaml title="params.yaml"
input: samplesheet.csv
outdir: results/
fasta: /path/to/genome.fa
gtf: /path/to/annotation.gtf
circrna_db: /path/to/circ_db.bed
circnick_species: mouse
```

---

## Pipeline parameters

### Reference files

| Parameter         | Description                                                           | Required |
| ----------------- | --------------------------------------------------------------------- | -------- |
| `--fasta`         | Reference genome FASTA file                                           | Yes      |
| `--gtf`           | Gene annotation GTF file                                              | Yes      |
| `--circrna_db`    | circRNA database BED file (required for isoCirc and CIRI-long)        | Conditional |
| `--genome_index_dir` | Directory to cache BWA genome index. Defaults to `<outdir>/genome_index` | No |

### Tool selection

By default all four tools are enabled, except that `--wet_lab circnick` runs CircNick-LRS alone and `--wet_lab isocirc` turns CIRI-long off (see below). Any `--run_*` flag overrides this. Disable individual tools with `false`:

| Parameter       | Description               | Default |
| --------------- | ------------------------- | ------- |
| `--run_isocirc`  | Run isoCirc               | `true` (`false` if `--wet_lab circnick`) |
| `--run_circfl`   | Run CircFL-seq            | `true` (`false` if `--wet_lab circnick`) |
| `--run_cirilong` | Run CIRI-long             | `true` (`false` if `--wet_lab circnick` or `isocirc`) |
| `--run_circnick` | Run circnick-lrs          | `true`  |

| Parameter    | Description | Default |
| ------------ | ----------- | ------- |
| `--wet_lab`  | Wet-lab protocol preset: `ciri_long`, `circfl`, `isocirc`, `circnick`. On real isoCirc-protocol data, CIRI-long has hung for 24+ days without finishing; `--wet_lab isocirc` defaults it off (`--run_cirilong true` opts back in). `--wet_lab circnick` is the only preset that changes tool selection beyond that. | `null` |

### circnick-lrs options

| Parameter                  | Description                                                                                     | Required when           |
| -------------------------- | ----------------------------------------------------------------------------------------------- | ----------------------- |
| `--circnick_species`       | Species for circnick-lrs built-in reference: `mouse` or `human`                                | `--run_circnick true`   |
| `--circnick_liftover_chain`| UCSC `.chain` file to lift circnick coordinates to the current genome build (optional)          | Never (optional)        |

> [!NOTE]
> circnick-lrs uses built-in mm10 (mouse) or hg19 (human) references internally. If your analysis uses a different genome build, provide `--circnick_liftover_chain` to convert coordinates.

### Tool-native options

Pass extra flags straight through to isoCirc's or CircFL-seq's own CLI, on top of nanocirc's own defaults (empty by default: each tool's own defaults apply). CIRI-long and circnick-lrs have no user-facing detection-tuning flags at all, so there is no equivalent flag for either.

| Parameter        | Description                                                              | Default |
| ----------------- | --------------------------------------------------------------------------- | ------- |
| `--isocirc_args` | Extra args to the `isocirc` command. Cannot include `-t`/`--threads`/`--bedtools`/`--minimap2`, which nanocirc already manages. | `''` |
| `--circfl_args`  | Extra args for circfull's `RG` step. Only `-u` (disable fusion-circRNA detection) and `-m <rmsk>` (repeat-region .bed path must be provided with it) are accepted; every other circfull flag is a hardcoded I/O path nanocirc uses to chain `RG`/`DNSC`/`cRG`/`mRG` together, so anything else is rejected at startup. | `''` |

```bash
nextflow run aerusakovich/nanocirc_pipeline \
    --input samplesheet.csv \
    --isocirc_args '--min-copy 3.0' \
    --circfl_args '-u' \
    ...
```

### Merge options

These options control how results from multiple tools are merged and scored:

| Parameter                   | Description                                                              | Default |
| --------------------------- | ------------------------------------------------------------------------ | ------- |
| `--circrna_bsj_tolerance`   | BSJ coordinate tolerance in bp for relaxed BSJ grouping during merge    | `5`     |
| `--circrna_isoform_overlap` | Minimum reciprocal spliced-length overlap for isoform confidence scoring | `0.95`  |
| `--circrna_trusted_tools`   | Comma-separated tools whose Low-confidence single-tool calls `balanced_recall`'s `trusted_only` filter still rescues | `cirilong,isocirc,circfl` |

Merging is only performed when **two or more** detection tools are active.

### Annotation options

| Parameter            | Description                                          | Default |
| --------------------- | ----------------------------------------------------- | ------- |
| `--skip_annotation`   | Skip GFF/FASTA/annotated TSV generation               | `false` |

### Cross-run merge

When the same sample (or biological condition) is sequenced across multiple runs, results can be merged across runs to increase sensitivity and confidence.

To enable cross-run merge, add a `group` column to the samplesheet. Samples sharing the same group name are merged together:

```csv
sample,fastq,group
run1,/data/run1.fq.gz,condition_A
run2,/data/run2.fq.gz,condition_A
run3,/data/run3.fq.gz,condition_A
ctrl1,/data/ctrl1.fq.gz,condition_B
ctrl2,/data/ctrl2.fq.gz,condition_B
```

Then set `--run_crossrun_merge true`:

```bash
nextflow run aerusakovich/nanocirc_pipeline \
    --input samplesheet.csv \
    --run_crossrun_merge true \
    ...
```

| Parameter              | Description                                           | Default |
| ---------------------- | ----------------------------------------------------- | ------- |
| `--run_crossrun_merge` | Enable cross-run merging using the `group` column     | `false` |
| `--crossrun_min_tool_agreement` | Min tool agreement needed to keep a single-run structure call | `2` |

> [!NOTE]
> The `group` column is optional. Samples without a group are processed per-sample only and are not included in any cross-run merge output.

### Quantification

Set `--run_quantify true` to enable circRNA read quantification.

Without `--run_crossrun_merge`, each sample is quantified against its own discovery catalog. With `--run_crossrun_merge`, every sample sharing a `group` is quantified against that group's crossrun discovery catalog instead.

```bash
nextflow run aerusakovich/nanocirc_pipeline \
    --input samplesheet.csv \
    --run_quantify true \
    ...
```

| Parameter                       | Description                                                       | Default |
| -------------------------------- | ------------------------------------------------------------------ | ------- |
| `--run_quantify`                 | Enable quantification                                              | `false` |
| `--quant_chunk_size`             | Loci per reference chunk for tier1                                 | `100`   |
| `--quant_chunk_seed`             | Fixed shuffle seed for chunk assignment                             | `42`    |
| `--quant_locus_dedup_tolerance`  | bp tolerance (both start AND end) for near-duplicate locus dedup   | `10`    |
| `--quant_min_old_tool_count`     | Min independent old-tool-count to flag a locus for tier2 rescue    | `50`    |
| `--circrna_confident_min_reads`  | Max quantified read count for the discovery/balanced_recall confidence filter (see above) | `2` |

`--run_quantify true` also builds a DESeq2-ready count matrix per confidence tier automatically; no additional flags. See [docs/output.md](output.md#deseq2-count-matrix) for the output files.

### QC options

| Parameter        | Description                      | Default |
| ---------------- | -------------------------------- | ------- |
| `--skip_qc`      | Skip all QC steps                | `false` |
| `--skip_fastqc`  | Skip FastQC                      | `false` |
| `--skip_nanoplot`| Skip NanoPlot                    | `false` |
| `--skip_multiqc` | Skip MultiQC report generation   | `false` |

---

## Profiles

Use `-profile` to configure the execution environment. Multiple profiles can be combined, e.g. `-profile singularity,genouest`.

| Profile       | Description                                      |
| ------------- | ------------------------------------------------ |
| `docker`      | Run with Docker containers                       |
| `singularity` | Run with Singularity containers                  |
| `apptainer`   | Run with Apptainer containers                    |
| `test`        | Minimal test run with test data          |

> [!IMPORTANT]
> This pipeline requires containers (Docker, Singularity, or Apptainer). The four circRNA detection tools are only available as container images and **conda is not supported**.

---

## Core Nextflow arguments

> [!NOTE]
> These use a single hyphen (`-`), unlike pipeline parameters which use double hyphen (`--`).

### `-resume`

Restart a pipeline reusing cached results where inputs are unchanged.

### `-work-dir` / `-w`

Directory for Nextflow working files. Defaults to `./work`. On HPC systems it is recommended to set this to a fast scratch filesystem.

### `-c`

Provide a custom Nextflow config file for tuning resource requirements or infrastructure settings.

---

## Resource requirements

Default resource labels used by the pipeline:

| Label              | CPUs | Memory  | Time   |
| ------------------ | ---- | ------- | ------ |
| `process_single`   | 1    | 6 GB    | 4 h    |
| `process_low`      | 2    | 12 GB   | 4 h    |
| `process_medium`   | 6    | 36 GB   | 8 h    |
| `process_high`     | 12   | 72 GB   | 16 h   |
| `process_long`     | 1 (default) | 6 GB (default) | 20 h |
| `process_high_memory` | 6 | 200 GB  | 4 h (default) |

`process_long` only overrides time; unset fields fall back to the pipeline-wide default (1 CPU, 6 GB, 4 h). It isn't used by any process by default but is available for your own config. To override resources for a specific process, add to your config:

```groovy
process {
    withName: 'ISOCIRC' {
        cpus   = 16
        memory = '100.GB'
        time   = '24.h'
    }
}
```

---

## Running on a cluster

By default Nextflow uses its `local` executor: every process runs as a
subprocess on whichever machine you launch `nextflow run` from.

**Check for an existing institution profile first.** If your cluster is listed at
[nf-core/configs](https://github.com/nf-core/configs/tree/master/conf), just add it
to `-profile` (e.g. `-profile apptainer,genouest`) and skip the config file below.
It already sets the right executor, queue, and resource limits for that cluster.

If there is no institution profile, tell Nextflow to submit each process as its own job
via a custom config (`-c`):

```groovy title="cluster.config"
process {
    executor = 'slurm'          // or 'sge', 'pbspro', 'lsf', ... see the Nextflow executor docs
    queue    = 'your_partition' // find with `sinfo -s`
}
```

Nextflow then submits every process (FastQC, each detection tool, quantification, ...)
as its own scheduler job, sized from the [resource labels](#resource-requirements)
above. The `nextflow run` command itself becomes a lightweight, long-lived orchestrator. See [Running in the background](#running-in-the-background)
below for keeping it alive.

> [!WARNING]
> Without either of the above, every process shares one job's walltime and CPU
> pool (Nextflow's `local` executor default). On real data, slower tools like
> `CIRI_LONG` and `CIRCNICK_LRS` can still be running when that job's time limit
> hits, killing the whole run. Symptoms: the session dies close to a round
> walltime boundary, a sibling process shows `exit: 143` (SIGTERM), and
> Nextflow may log `process hasn't exited` right before aborting. That is the
> outer job's walltime expiring, not a pipeline bug.

---

## Running in the background

Use `screen`, `tmux`, or the Nextflow `-bg` flag to detach the run from your terminal:

```bash
nextflow run aerusakovich/nanocirc_pipeline ... -bg
```

Alternatively, on HPC systems, submit the Nextflow head job itself to the scheduler:

```bash
sbatch --wrap="nextflow run aerusakovich/nanocirc_pipeline ..."
```

---

## Troubleshooting

### Pipeline exits with "At least one tool must be active"

All four tool flags are `true` by default. This error only appears if you explicitly set all of them to `false`. Enable at least one tool.

### Pipeline exits with "CircRNA analysis requires '--circrna_db'"

isoCirc and CIRI-long require a circRNA database BED file. Either provide `--circrna_db` or disable those tools with `--run_isocirc false --run_cirilong false`.

### circnick-lrs exits with coordinate warnings

If many circRNAs report exons outside their BSJ boundaries, consider providing `--circnick_liftover_chain` to lift coordinates to the current genome build.

### Memory errors on genome indexing

If PREPARE_GENOME fails due to memory, override its resources in your config:

```groovy
process {
    withName: 'PREPARE_GENOME' {
        memory = '32.GB'
    }
}
```

### Memory errors from CIRI-long or CircFL-seq

Both tools multi-thread by splitting work across cores, and each core's share adds to peak memory, so more CPUs does not always mean faster-and-safer - it can push memory usage past what's available. If one of these runs out of memory, try reducing its `cpus` rather than only raising `memory`:

```groovy
process {
    withName: 'CIRI_LONG' {
        cpus = 4
    }
    withName: 'CIRCFL_SEQ' {
        cpus = 4
    }
}
```

For output file formats see [docs/output.md](output.md). For the merge/scoring algorithms see [docs/methods.md](methods.md).
