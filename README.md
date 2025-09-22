# ATAC-seq Analysis Pipeline

This repository contains a **generalized ATAC-seq processing pipeline**.  
It is designed to take raw FASTQ files and produce **cleaned, aligned, deduplicated BAM files** ready for downstream peak calling and analysis.

---

## 🚀 Workflow Overview
1. **Quality Control (FastQC)** – Assess raw FASTQ quality.
2. **Adapter Trimming (fastp)** – Remove adapters and low-quality reads.
3. **QC after trimming** – Verify trimming effectiveness.
4. **Alignment (Bowtie2)** – Map reads to a reference genome.
5. **Sorting & Indexing (Samtools)** – Prepare BAM for downstream steps.
6. **Mitochondrial Read Removal** – Filter out chrM reads.
7. **Deduplication (Picard)** – Remove PCR duplicates.

---

## 🔧 Requirements
Install the following tools before running the pipeline:
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [fastp](https://github.com/OpenGene/fastp)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/)
- [Samtools](http://www.htslib.org/)
- [Picard](https://broadinstitute.github.io/picard/)

---

## ⚙️ Usage
```bash
bash ATAC_seq.sh <INPUT_DIR> <OUTPUT_DIR> <BOWTIE2_INDEX>
