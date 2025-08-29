

mkdir -p fastqc && fastqc -t $(nproc) G2025_9/*_R[12]_001.fastq.gz -d fastqc -o fastqc


mkdir -p trimmed && \
ls G2025_9/*_R1_001.fastq.gz | \
xargs -n1 -P $(( $(nproc)/16 )) -I{} bash -c '
  s=$(basename {} _R1_001.fastq.gz)
  fastp -w 16 \
        -i "{}" \
        -I "G2025_9/${s}_R2_001.fastq.gz" \
        -o "trimmed/${s}_R1.trimmed.fastq.gz" \
        -O "trimmed/${s}_R2.trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        -j "trimmed/${s}.fastp.json" \
        -h "trimmed/${s}.fastp.html"
'

mkdir -p fastqc_trimmed && fastqc -t $(nproc) trimmed/*_R[12].trimmed.fastq.gz -d fastqc -o fastqc_trimmed


mkdir -p aligned && \
for f in trimmed/*_R1.trimmed.fastq.gz; do \
  s=$(basename "$f" _R1.trimmed.fastq.gz); \
  bowtie2 -p $(nproc) \
    --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 \
    -x /data/shared_genomics_data/ref_genomes/human/hg38/bowtie2/hg38 \
    -1 "$f" \
    -2 "trimmed/${s}_R2.trimmed.fastq.gz" \
  | samtools view -@ $(nproc) -bS - > "aligned/${s}.bam"; \
done



##sort
mkdir -p sorted && \
for f in aligned/*.bam; do \
  s=$(basename "$f" .bam); \
  samtools sort -@ $(nproc) -o sorted/${s}_sorted.bam "$f" && \
  samtools index -@ $(nproc) sorted/${s}_sorted.bam; \
done



##check chrM
mkdir -p idxstats && \
ls sorted/*_sorted.bam | xargs -P $(nproc) -I {} bash -c '
  s=$(basename {} _sorted.bam)
  samtools idxstats {} > idxstats/${s}_sorted.idxstats
  grep "chrM" idxstats/${s}_sorted.idxstats > idxstats/${s}_sorted.chrM.txt
'


##remove chrM
mkdir -p rmChrM && \
for f in sorted/*_sorted.bam; do \
  s=$(basename "$f" _sorted.bam); \
  echo "=== Processing $s ==="; \
  samtools view -h "$f" | grep -v chrM | \
    samtools sort -@ $(nproc) -O bam -T rmChrM/${s}_tmp -o rmChrM/${s}.rmChrM.bam && \
  echo "Created rmChrM/${s}.rmChrM.bam"; \
done




##dedup

mkdir -p dedup && \
for f in sorted/*_sorted.bam; do \
  s=$(basename "$f" _sorted.bam); \
  nohup picard AddOrReplaceReadGroups \
    I="$f" \
    O="dedup/${s}_sorted_rg.bam" \
    RGID="$s" \
    RGLB="lib_${s}" \
    RGPL="illumina" \
    RGPU="unit_${s}" \
    RGSM="${s}_sample" \
  > "dedup/${s}_redup.log" 2>&1 & \
done










