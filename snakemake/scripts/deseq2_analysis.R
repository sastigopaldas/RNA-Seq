#!/usr/bin/env Rscript
# DESeq2 Differential Expression Analysis Script
# This script performs differential expression analysis using DESeq2

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(EnhancedVolcano)
    library(ggrepel)
})

# ---------------------------
# Get input/output parameters
# ---------------------------
counts_file <- snakemake@input[["counts"]]
samples_file <- snakemake@input[["samples"]]

results_file <- snakemake@output[["results"]]
normalized_counts_file <- snakemake@output[["normalized_counts"]]
rlog_counts_file <- snakemake@output[["rlog_counts"]]
pca_plot_file <- snakemake@output[["pca_plot"]]
ma_plot_file <- snakemake@output[["ma_plot"]]
volcano_plot_file <- snakemake@output[["volcano_plot"]]
heatmap_file <- snakemake@output[["heatmap"]]

contrast <- snakemake@params[["contrast"]]     # e.g., c("condition", "treated", "control")
alpha <- snakemake@params[["alpha"]]           # e.g., 0.05
lfc_threshold <- snakemake@params[["lfc_threshold"]]  # e.g., 1
cat("Loading data...\n")

# ---------------------------
# Load and prepare data
# ---------------------------
sample_info <- read_tsv(samples_file, show_col_types = FALSE)
sample_info <- sample_info[!duplicated(sample_info$sample), ]
rownames(sample_info) <- sample_info$sample

count_data <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Remove first 5 columns (featureCounts metadata)
count_data <- count_data[, -(1:5)]

# Clean sample names
colnames(count_data) <- gsub(".*/", "", colnames(count_data))
colnames(count_data) <- gsub("\\.sorted\\.bam$", "", colnames(count_data))

# Ensure matching samples
common_samples <- intersect(colnames(count_data), rownames(sample_info))
count_data <- count_data[, common_samples]
sample_info <- sample_info[common_samples, ]

cat(paste("Found", ncol(count_data), "samples and", nrow(count_data), "genes\n"))

# Explicitly set factors
sample_info$condition <- factor(sample_info$condition, levels = c("control", "treated"))
sample_info$batch <- factor(sample_info$batch)  # optional if adjusting for batch

# ---------------------------
# Create DESeq dataset
# ---------------------------
dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_info,
    design = ~ batch + condition  # include batch or just ~ condition
)




# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(paste("After filtering:", nrow(dds), "genes retained\n"))

# ---------------------------
# Run DESeq2
# ---------------------------
cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)

res <- results(dds,
               contrast = c("condition","treated","control"),
               alpha = alpha,
               lfcThreshold = lfc_threshold)

sig_genes <- sum(res$padj < alpha & abs(res$log2FoldChange) > lfc_threshold, na.rm = TRUE)
cat(paste("Found", sig_genes, "significantly differentially expressed genes\n"))

# ---------------------------
# Save results
# ---------------------------
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write_csv(res_df, results_file)

# Normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$gene_id <- rownames(normalized_counts_df)
write_csv(normalized_counts_df, normalized_counts_file)

# ---------------------------
# rlog transformation
# ---------------------------
cat("Performing rlog transformation...\n")
rld <- rlog(dds, blind = FALSE)
rlog_counts_df <- as.data.frame(assay(rld))
rlog_counts_df$gene_id <- rownames(rlog_counts_df)
write_csv(rlog_counts_df, rlog_counts_file)

# ---------------------------
# PCA plot
# ---------------------------
cat("Generating PCA plot...\n")
pca_data <- plotPCA(rld, intgroup = contrast[1], returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = !!sym(contrast[1]))) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = name), size = 3) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    theme_minimal() +
    ggtitle("Principal Component Analysis")
ggsave(pca_plot_file, p_pca, width = 8, height = 6, dpi = 300)

# ---------------------------
# MA plot
# ---------------------------
cat("Generating MA plot...\n")
png(ma_plot_file, width = 2400, height = 1800, res = 300)
plotMA(res, alpha = alpha, ylim = c(-5, 5))
title(main = paste("MA Plot -", contrast[2], "vs", contrast[3]))
dev.off()

# ---------------------------
# Volcano plot
# ---------------------------
cat("Generating volcano plot...\n")
volcano_plot <- EnhancedVolcano(
    res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "padj",
    title = paste("Volcano Plot -", contrast[2], "vs", contrast[3]),
    pCutoff = alpha,
    FCcutoff = lfc_threshold,
    pointSize = 2.0,
    labSize = 3.0,
    colAlpha = 1,
    legendPosition = "right",
    drawConnectors = TRUE
)
ggsave(volcano_plot_file, volcano_plot, width = 10, height = 8, dpi = 300)

# ---------------------------
# Sample correlation heatmap
# ---------------------------
cat("Generating sample correlation heatmap...\n")
sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- paste(rld$sample, rld[[contrast[1]]], sep = " - ")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(heatmap_file, width = 2400, height = 2400, res = 300)
pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    col = colors,
    main = "Sample-to-Sample Distances"
)
dev.off()

cat("✅ DESeq2 analysis completed successfully!\n")

