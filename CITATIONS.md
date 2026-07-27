# nanocirc: Citations

## nanocirc

If you use this pipeline, please cite:

> Rusakovich A, Derrien T*, Blum Y* (*co-last authors). nanocirc: a Nextflow pipeline for long-read circRNA detection and annotation using a consensus multi-tool approach and confidence scoring. Manuscript in preparation.

This pipeline's design and benchmark was motivated by our first benchmark study:

> Rusakovich A, Corre S, Cadieu E, Fraboulet RM, Le Bars V, Galibert MD, Derrien T, Blum Y. Benchmarking circRNA detection tools from long-read sequencing using a data-driven and flexible simulation framework. Peer Community Journal. 2026;6:e27. doi: 10.24072/pcjournal.699.

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

## Pipeline tools

### circRNA detection

- [isoCirc](https://github.com/Xinglab/isoCirc)

  > Xin R, Gao Y, Gao Y, Wang R, Kadash-Edmondson KE, Liu B, Wang Y, Lin L, Xing Y. isoCirc catalogs full-length circular RNA isoforms in human transcriptomes. Nat Commun. 2021 Jan 12;12(1):266. doi: 10.1038/s41467-020-20459-8. PubMed PMID: 33436586; PubMed Central PMCID: PMC7803974.

- [CircFL-seq](https://github.com/yangence/circfull)

  > Liu Z, Tao C, Li S, Du M, Bai Y, Hu X, Li Y, Chen J, Yang E. circFL-seq reveals full-length circular RNAs with rolling circular reverse transcription and nanopore sequencing. eLife. 2021 Nov 16;10:e69457. doi: 10.7554/eLife.69457. PubMed PMID: 34783652; PubMed Central PMCID: PMC8594938.

- [CIRI-long](https://github.com/bioinfo-biols/CIRI-long)

  > Zhang J, Hou L, Zuo Z, Ji P, Zhang X, Xue Y, Zhao F. Comprehensive profiling of circular RNAs with nanopore sequencing and CIRI-long. Nat Biotechnol. 2021 Jul;39(7):836-845. doi: 10.1038/s41587-021-00842-6. PubMed PMID: 33782612.

- [circnick-lrs](https://github.com/dzhang32/circnick)

  > Rahimi K, Venø MT, Dupont DM, Kjems J. Nanopore sequencing of brain-derived full-length circRNAs reveals circRNA-specific exon usage, intron retention and microexons. Nat Commun. 2021 Aug 10;12(1):4825. doi: 10.1038/s41467-021-24975-z. PubMed PMID: 34376670; PubMed Central PMCID: PMC8355219.

### Merging and comparison

- [bedtools](https://pubmed.ncbi.nlm.nih.gov/20110278/)

  > Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033. PubMed PMID: 20110278; PubMed Central PMCID: PMC2832824.

- [pybedtools](https://pubmed.ncbi.nlm.nih.gov/21949271/)

  > Dale RK, Pedersen BS, Quinlan AR. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. Bioinformatics. 2011 Dec 15;27(24):3423-4. doi: 10.1093/bioinformatics/btr539. PubMed PMID: 21949271; PubMed Central PMCID: PMC3232374.

### Coordinate liftover

- [UCSC liftOver / Genome Browser](https://pubmed.ncbi.nlm.nih.gov/12045153/)

  > Kent WJ, Sugnet CW, Furey TS, Roskin KM, Pringle TH, Zahler AM, Haussler D. The human genome browser at UCSC. Genome Res. 2002 Jun;12(6):996-1006. doi: 10.1101/gr.229102. PubMed PMID: 12045153; PubMed Central PMCID: PMC186604.

### Annotation

- [samtools](https://pubmed.ncbi.nlm.nih.gov/33590861/)

  > Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PubMed PMID: 33590861; PubMed Central PMCID: PMC7931819.

- [GffRead and GffCompare](https://pubmed.ncbi.nlm.nih.gov/32489650/)

  > Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. F1000Res. 2020 Apr 28;9:ISCB Comm J-304. doi: 10.12688/f1000research.23297.2. PubMed PMID: 32489650; PubMed Central PMCID: PMC7222033.

- [AGAT](https://www.doi.org/10.5281/zenodo.3552717)

  > Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. (Version v1.0.0). Zenodo. doi: 10.5281/zenodo.3552717.

### Quantification

- [minimap2](https://pubmed.ncbi.nlm.nih.gov/29750242/)

  > Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191. PubMed PMID: 29750242; PubMed Central PMCID: PMC6137996.

- [BLAT](https://pubmed.ncbi.nlm.nih.gov/11932250/)

  > Kent WJ. BLAT--the BLAST-like alignment tool. Genome Res. 2002 Apr;12(4):656-64. doi: 10.1101/gr.229202. PubMed PMID: 11932250; PubMed Central PMCID: PMC187518.

- [BWA](https://pubmed.ncbi.nlm.nih.gov/19451168/)

  > Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics. 2009 Jul 15;25(14):1754-60. doi: 10.1093/bioinformatics/btp324. PubMed PMID: 19451168; PubMed Central PMCID: PMC2705234.

### Quality control

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

  > Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].

- [NanoPlot](https://pubmed.ncbi.nlm.nih.gov/29547981/)

  > De Coster W, D'Hert S, Schultz DT, Cruts M, Van Broeckhoven C. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. 2018 Aug 1;34(15):2666-2669. doi: 10.1093/bioinformatics/bty149. PubMed PMID: 29547981; PubMed Central PMCID: PMC6061794.

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

## Software packaging/containerisation tools

- [containers_circRNA_tools](https://github.com/aerusakovich/containers_circRNA_tools) - this pipeline's own patched container builds for isoCirc, CircFL-seq, CIRI-long, circnick-lrs and nanocirc.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.
