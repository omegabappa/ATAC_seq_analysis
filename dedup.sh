#!/usr/bin/env bash
set -euo pipefail
trap 'echo "❌ Error on line $LINENO"; exit 1' ERR

threads=128
outdir="dedup"

echo "Starting ReadGroup+MarkDuplicates on all samples with $threads threads…"
mkdir -p "$outdir"

for bam in rmChrM/*_rmChrM.bam; do
  sample="$(basename "$bam" _rmChrM.bam)"
  echo
  echo "⏳ Processing $sample …"

  # 1) Add/Replace Read Groups
  echo "   → AddOrReplaceReadGroups"
  java -XX:ActiveProcessorCount=$threads -jar /home/bappa.ghosh/bin/picard.jar AddOrReplaceReadGroups \
    I="$bam" \
    O="$outdir/${sample}_rg.bam" \
    RGID="$sample" \
    RGLB="${sample}_lib" \
    RGPL="illumina" \
    RGPU="unit1" \
    RGSM="$sample"

  echo "     ✔ RG added: $outdir/${sample}_rg.bam"

  # 2) MarkDuplicates (and remove them)
  echo "   → MarkDuplicates"
  java -XX:ActiveProcessorCount=$threads -jar /home/bappa.ghosh/bin/picard.jar MarkDuplicates \
    I="$outdir/${sample}_rg.bam" \
    O="$outdir/${sample}_dedup.bam" \
    M="$outdir/${sample}_dup_metrics.txt" \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true

  echo "     ✔ Duplicates removed: $outdir/${sample}_dedup.bam"
done

echo
echo "✅ All samples processed. Check '$outdir/' for outputs and metrics."

