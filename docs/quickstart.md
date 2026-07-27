# nanocirc: Quick Start Guide

A simple guide for your first real run. If you're comfortable with Nextflow already, [docs/usage.md](usage.md) is the faster reference.

## 1. Install what you need

- **Java 17 or newer** (Nextflow's own requirement; this pipeline is developed against Java 21). If you don't have Java, [SDKMAN](https://sdkman.io/) is a simple way to install one: `curl -s "https://get.sdkman.io" | bash`, then `sdk install java 21.0.3-tem`.
- **Nextflow >=25.04.0** (this pipeline's minimum, see `nextflow.config`). Install with `curl -s https://get.nextflow.io | bash`, or see the [official install guide](https://nf-co.re/docs/usage/installation).
- **A container runtime**: Singularity, Apptainer, or Docker. The four circRNA detection tools only come as containers, conda is not supported due to reproducibility reasons. On a shared HPC cluster, Singularity/Apptainer is usually already installed system-wide, check with your cluster admin or `module avail` if unsure.

## 2. Files you need

| File | What it is | Required when |
| ---- | ---------- | -------------- |
| FASTQ | Your own nanopore long-read data, gzipped (`.fastq.gz`/`.fq.gz`) | Always |
| Genome FASTA | Reference genome matching your reads | Always |
| GTF | Gene annotation matching the same genome build as the FASTA | Always |
| circRNA database BED | Known circRNA loci to check candidates against | `--run_isocirc` or `--run_cirilong` (both default `true`) |
| Liftover chain | Converts circnick-lrs' output coordinates to your genome build | If your genome version isn't circnick-lrs' built-in hg19 (human) or mm10 (mouse) |

For the genome FASTA and GTF, [GENCODE](https://www.gencodegenes.org/human/) (human) and the [mouse GENCODE page](https://www.gencodegenes.org/mouse/) provide both files together, matched to the same build, which avoids a common source of coordinate mismatches. [UCSC's goldenPath](https://hgdownload.soe.ucsc.edu/goldenPath/) is another standard source.

## 3. circRNA database BED file

isoCirc and CIRI-long both check candidate back-splice junctions against a reference circRNA database. It must be tab-separated BED4, no header:

```
chrom  start  end  strand
```

We provide two ready-to-use circRNA databases based on CircBase+CircAtlas :

- `assets/circrna_db/circrna_human_db.bed` (639,211 loci, ~16 MB)
- `assets/circrna_db/circrna_mouse_db.bed` (344,719 loci, ~15 MB)

You can point your `--circrna_db` straight at whichever matches your species:

```bash
--circrna_db assets/circrna_db/circrna_human_db.bed
```

For chr21-only sanity-check or as a formatting template, there's also a smaller subset: `assets/test_data/human_chr21/chr21_circrna_db.bed`. For other species or an updated database, [circBase](http://www.circbase.org/) and circAtlas are common public sources for circRNA coordinates; export/reformat to BED4 before use.

## 4. Liftover chain 

circnick-lrs only comes with a built-in hg19 (human) or mm10 (mouse) reference. If your own FASTA/GTF use a different build (e.g. hg38), provide `--circnick_liftover_chain` so circnick-lrs' output gets lifted to match everything else.

Standard UCSC chain files are available under `https://hgdownload.soe.ucsc.edu/goldenPath/<build>/liftOver/`, e.g. [hg19 chains](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) for human, [mm10 chains](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/) for mouse.

We demonstrate a working example: `assets/test_data/human_chr21/hg19ToHg38.over.chain.gz`.

## 5. Samplesheet

CSV, one row per FASTQ file:

```csv
sample,fastq
SAMPLE1,/path/to/sample1.fastq.gz
SAMPLE2,/path/to/sample2.fastq.gz
```

Add a `group` column if you plan to use `--run_crossrun_merge` (merges circRNA calls across sequencing runs of the same biological sample/condition):

```csv
sample,fastq,group
run1,/data/run1.fastq.gz,condition_A
run2,/data/run2.fastq.gz,condition_A
```

## 6. Config: local machine vs cluster

**Local machine** (laptop/workstation): no config file needed, just pick a container profile:

```bash
nextflow run aerusakovich/nanocirc_pipeline -profile singularity --input samplesheet.csv --outdir results/ ...
```

**HPC cluster** (e.g. SLURM): add a config file telling Nextflow to submit each process as its own job, and where your container cache is located:

```groovy
// cluster.config
process {
    executor = 'slurm'
}
singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = '/path/to/shared/singularity_cache'
}
```

Then add `-c cluster.config` to the command:

```bash
nextflow run aerusakovich/nanocirc_pipeline -profile singularity -c cluster.config --input samplesheet.csv --outdir results/ ...
```

See [docs/usage.md](usage.md#resource-requirements) to override CPU/memory/time for specific processes.

## 7. Run file example

Put the whole command in a script so the run is reproducible and easy to resubmit with `-resume`:

```bash
#!/bin/bash
#SBATCH --job-name=nanocirc_run
#SBATCH --output=nanocirc_%j.log
#SBATCH --time=120:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

nextflow run aerusakovich/nanocirc_pipeline \
    -profile singularity \
    -c cluster.config \
    --input             samplesheet.csv \
    --outdir            results/ \
    --fasta             genome.fa \
    --gtf               annotation.gtf \
    --circrna_db        circrna_db.bed \
    --circnick_species  human \
    -resume
```

The `#SBATCH` lines above are only for the Nextflow *head* job, a lightweight coordinator that itself submits each pipeline step as its own separate job. If you don't use SLURM, drop them and just run the script in a `screen` or `tmux` session so it survives you disconnecting.

## 8. Test your setup first

Before trying real data, confirm everything works with the provided test dataset, no samplesheet or reference files needed:

```bash
nextflow run aerusakovich/nanocirc_pipeline -profile test,singularity --outdir test_results/
```

If it finishes with `Pipeline completed successfully`, your installation is good to go.

---

For full parameter reference, see [docs/usage.md](usage.md). For what the pipeline produces, see [docs/output.md](output.md).