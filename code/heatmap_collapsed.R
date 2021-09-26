setwd('Desktop/PUA_LAB/neil')

info_data <- read.csv("neil/deseq2_collapsed/HFD_ff_v_HFD_ko.csv",header=T,row.names=1)

info_data$log2FoldChange <- info_data$log2FoldChange * -1

norm_counts <- read.csv("neil/deseq2_collapsed/norm_counts.csv",header=T,row.names=1)
total_reads <- colSums(norm_counts)
norm_counts <- norm_counts/total_reads*1000000
norm_counts <- log2(norm_counts)

upreg <- subset(info_data, baseMean > 100 & log2FoldChange > 0.7 & padj < 0.05)
norm_counts$names <- rownames(norm_counts)
upreg$names <- rownames(upreg)
up_norm <- subset(norm_counts, norm_counts$names %in% upreg$names, select = c(X1,X2,X3,X4,X5,X6,X7,X8))
up_norm <- as.matrix(up_norm)

downreg <- subset(info_data, baseMean > 100 & log2FoldChange < -0.7 & padj < 0.05)
downreg$names <- rownames(downreg)
down_norm <- subset(norm_counts, norm_counts$names %in% downreg$names, select = c(X1,X2,X3,X4,X5,X6, X7,X8))
down_norm <- as.matrix(down_norm)

library(pheatmap)
pheatmap(up_norm, scale = "none", cellwidth = 20, cellheight = 15, cluster_rows = TRUE, cluster_cols = FALSE, labels_col = c("HFD_ff_1", "HFD_ff_2", "HFD_ff_3", "HFD_ff_4", "HFD_ko_5", "HFD_ko_6", "HFD_ko_7", "HFD_ko_8"))
pheatmap(down_norm, scale = "none", cellwidth = 20, cellheight = 10, cluster_rows = TRUE, cluster_cols = FALSE, labels_col = c("HFD_ff_1", "HFD_ff_2", "HFD_ff_3", "HFD_ff_4", "HFD_ko_5", "HFD_ko_6", "HFD_ko_7", "HFD_ko_8"))

