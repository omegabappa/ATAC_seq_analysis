#!/bin/bash
# ==========================================================
# Generalized ATAC-seq Processing Pipeline
# Author: Your Name
# Description:
#   This script performs a standard ATAC-seq data processing 
#   workflow starting from raw FASTQ files.
#   Steps include: QC, trimming, alignment, filtering, deduplication.
#
# Requirements:
#   - fastqc
#   - fastp
#   - bowtie2
#   - samtools
#   - picard
#
# Usage:
#   bash ATAC_seq.sh <INPUT_DIR> <OUTPUT_DIR> <BOWTIE2_INDEX>
# Example:
#   bash ATAC_seq.sh ./raw_data ./results /path/to/hg38/index
# ==========================================================

set -euo pipefail

# ---------------------------
# Input arguments
# ---------------------------
INPUT_DIR=$1        # Directory containing raw FASTQ files
OUTPUT_DIR=$2       # Directory where results will be stored
BOWTIE2_INDEX=$3    # Path to Bowtie2 genome index

THREADS=$(nproc)    # Use all available cores

# ---------------------------
# Step 1: Initial QC (FastQC)
# ---------------------------
mkdir -p ${OUTPUT_DIR}/fastqc
fastqc -t ${THREADS} ${INPUT_DIR}/*_R[12]_001.fastq.gz -o ${OUTPUT_DIR}/fastqc

# ---------------------------
# Step 2: Adapter trimming (fastp)
# ---------------------------
mkdir -p ${OUTPUT_DIR}/trimmed
for f in ${INPUT_DIR}/*_R1_001.fastq.gz; do
  sample=$(basename $f _R1_001.fastq.gz)
  fastp -w 16 \
        -i ${INPUT_DIR}/${sample}_R1_001.fastq.gz \
        -I ${INPUT_DIR}/${sample}_R2_001.fastq.gz \
        -o ${OUTPUT_DIR}/trimmed/${sample}_R1.trimmed.fastq.gz \
        -O ${OUTPUT_DIR}/trimmed/${sample}_R2.trimmed.fastq.gz \
        --detect_adapter_for_pe \
        -j ${OUTPUT_DIR}/trimmed/${sample}.fastp.json \
        -h ${OUTPUT_DIR}/trimmed/${sample}.fastp.html
done

# ---------------------------
# Step 3: QC after trimming
# ---------------------------
mkdir -p ${OUTPUT_DIR}/fastqc_trimmed
fastqc -t ${THREADS} ${OUTPUT_DIR}/trimmed/*_R[12].trimmed.fastq.gz -o ${OUTPUT_DIR}/fastqc_trimmed

# ---------------------------
# Step 4: Alignment (Bowtie2)
# ---------------------------
mkdir -p ${OUTPUT_DIR}/aligned
for f in ${OUTPUT_DIR}/trimmed/*_R1.trimmed.fastq.gz; do
  sample=$(basename $f _R1.trimmed.fastq.gz)
  bowtie2 -p ${THREADS} \
    --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 \
    -x ${BOWTIE2_INDEX} \
    -1 $f \
    -2 ${OUTPUT_DIR}/trimmed/${sample}_R2.trimmed.fastq.gz \
  | samtools view -@ ${THREADS} -bS - > ${OUTPUT_DIR}/aligned/${sample}.bam
done

# ---------------------------
# Step 5: Sort & index BAM
# ---------------------------
mkdir -p ${OUTPUT_DIR}/sorted
for f in ${OUTPUT_DIR}/aligned/*.bam; do
  sample=$(basename $f .bam)
  samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/sorted/${sample}_sorted.bam $f
  samtools index -@ ${THREADS} ${OUTPUT_DIR}/sorted/${sample}_sorted.bam
done

# ---------------------------
# Step 6: Mitochondrial reads check & removal
# ---------------------------
mkdir -p ${OUTPUT_DIR}/rmChrM
for f in ${OUTPUT_DIR}/sorted/*_sorted.bam; do
  sample=$(basename $f _sorted.bam)
  echo "Processing $sample..."
  samtools idxstats $f > ${OUTPUT_DIR}/rmChrM/${sample}.idxstats
  grep "chrM" ${OUTPUT_DIR}/rmChrM/${sample}.idxstats > ${OUTPUT_DIR}/rmChrM/${sample}.chrM.txt
  samtools view -h $f | grep -v chrM | \
    samtools sort -@ ${THREADS} -O bam -T ${OUTPUT_DIR}/rmChrM/${sample}_tmp -o ${OUTPUT_DIR}/rmChrM/${sample}.rmChrM.bam
done

# ---------------------------
# Step 7: Deduplication (Picard)
# ---------------------------
mkdir -p ${OUTPUT_DIR}/dedup
for f in ${OUTPUT_DIR}/rmChrM/*.bam; do
  sample=$(basename $f .bam)
  picard MarkDuplicates \
    I=$f \
    O=${OUTPUT_DIR}/dedup/${sample}_dedup.bam \
    M=${OUTPUT_DIR}/dedup/${sample}_dedup_metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT
  samtools index ${OUTPUT_DIR}/dedup/${sample}_dedup.bam
done

echo "=== ATAC-seq pipeline completed successfully ==="
