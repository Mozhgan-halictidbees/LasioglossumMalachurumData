# RNA-seq Processing Pipeline (Bash)

This repository contains a complete bash-based RNA-seq processing pipeline for quality control, rRNA removal, trimming, alignment, and gene expression quantification.

---

## Pipeline Overview

The workflow includes the following steps:

1. FastQC (Raw reads quality control)
2. SortMeRNA (rRNA removal)
3. FastQC (Post rRNA filtering QC)
4. Trimmomatic (Adapter and quality trimming)
5. FastQC (Post trimming QC)
6. HISAT2 (Alignment to reference genome)
7. Samtools (SAM to BAM conversion)
8. HTSeq-count (Gene expression quantification)

---

## Input Data

Paired-end FASTQ files are required.

Supported naming formats:
- `*_1.fq.gz` / `*_2.fq.gz`
- `*_fwd.fq.gz` / `*_rev.fq.gz`

---

## Requirements

The following tools must be installed:

- FastQC
- SortMeRNA
- Trimmomatic
- HISAT2
- Samtools
- HTSeq-count
- Bash
- awk, grep, bc

Most tools can be installed using conda:

```bash
conda create -n rnaseq_env -c bioconda fastqc sortmerna trimmomatic hisat2 samtools htseq
conda activate rnaseq_env
