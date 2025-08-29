# 1) load library
library(DESeq2)

# 2) read in your counts file (skip the “# Program…” comment line)
fc <- read.delim("peak_counts.txt", header=TRUE, comment.char="#", stringsAsFactors=FALSE)

# 3) extract counts matrix and set rownames
cts <- as.matrix(fc[ , 7:ncol(fc)] )
rownames(cts) <- fc$Geneid

# 4) rename columns to match your replicates
colnames(cts) <- c(
  "RPE1_rep1", "RPE1_rep2",           # sample1 → RPE1 (2 reps)
  "KO1_rep1",  "KO1_rep2",  "KO1_rep3",# sample2 → KO1 (3 reps)
  "KO2_rep1",  "KO2_rep2",  "KO2_rep3" # sample3 → KO2 (3 reps)
)

# 5) build your sample‐info (colData) table
condition <- factor(rep(c("RPE1","KO1","KO2"), times=c(2,3,3)))
colData <- data.frame(condition=condition, row.names=colnames(cts))

# 6) create DESeq2 object and run
dds <- DESeqDataSetFromMatrix(countData=cts, colData=colData, design=~condition)
dds <- DESeq(dds)

# 7) extract contrasts
res_KO1_vs_RPE1 <- results(dds, contrast=c("condition","KO1","RPE1"))
res_KO2_vs_RPE1 <- results(dds, contrast=c("condition","KO2","RPE1"))

# 8) view top hits
head(res_KO1_vs_RPE1)
head(res_KO2_vs_RPE1)

# 9) write out full result tables
write.csv(as.data.frame(res_KO1_vs_RPE1),
          file="DESeq2_KO1_vs_RPE1_results.csv", row.names=TRUE)
write.csv(as.data.frame(res_KO2_vs_RPE1),
          file="DESeq2_KO2_vs_RPE1_results.csv", row.names=TRUE)


png("MA_KO1_vs_RPE1.png", width=800, height=600)
plotMA(res_KO1_vs_RPE1, main="KO1 vs RPE1", ylim=c(-5,5))
dev.off()






# Load libraries
library(ggplot2)

# Helper to make & save a volcano
make_volcano <- function(res, filename, title) {
  df <- as.data.frame(res)
  
  # 1) drop NAs in padj or log2FC
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  
  # 2) replace padj==0 with smallest nonzero padj to avoid Inf
  if(any(df$padj == 0, na.rm=TRUE)) {
    min_nonzero <- min(df$padj[df$padj > 0], na.rm=TRUE)
    df$padj[df$padj == 0] <- min_nonzero
  }
  
  # flag significance
  df$signif <- with(df, padj < 0.05 & abs(log2FoldChange) > 1)
  
  # build plot
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(title = title,
         x = "log2 fold‐change",
         y = "-log10(padj)")
  
  # save as JPEG @600 dpi
  ggsave(
    filename = filename,
    plot     = p,
    width    = 6, 
    height   = 5, 
    units    = "in",
    dpi      = 600,
    device   = "jpeg"
  )
}

# Volcano 1: KO1 vs RPE1
make_volcano(
  res      = res_KO1_vs_RPE1,
  filename = "volcano_KO1_vs_RPE1.jpg",
  title    = "Volcano: KO1 vs RPE1"
)

# Volcano 2: KO2 vs RPE1
make_volcano(
  res      = res_KO2_vs_RPE1,
  filename = "volcano_KO2_vs_RPE1.jpg",
  title    = "Volcano: KO2 vs RPE1"
)







#4. Heatmap of top DE peaks

# variance‑stabilize
vsd <- vst(dds, blind=FALSE)

# pick top 50 by padj
top50 <- rownames(head(res_KO1_vs_RPE1[order(res_KO1_vs_RPE1$padj), ], 50))
mat <- assay(vsd)[ top50, ]

# annotation for columns
anno_col <- as.data.frame(colData(dds)["condition", drop=FALSE])

pheatmap(mat,
         cluster_rows=TRUE, cluster_cols=TRUE,
         annotation_col=anno_col,
         show_rownames=FALSE,
         filename="heatmap_top50_KO1_vs_RPE1.png",
         width=6, height=8)






library(DESeq2); library(ggplot2)
library(pheatmap); library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db); library(clusterProfiler)


