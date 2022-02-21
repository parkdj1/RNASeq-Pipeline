# ###########################
#
#  heatmap_collapsed.R
#
#  Created on Wed May 6 2020
#  @author: parkdj1
#
#  create heatmaps of up/down-regulated genes from collapsed deseq and norm counts data
#
#  input: paired reads CSV file | normalized counts CSV file
#  output: upregulated genes heatmap | downregulated genes heatmap
#
# ###########################

# set working directory (AKA where are your files going to be/output to)
setwd('../')

# set path to CSV files containing data
deseq_file <- "PATH/TO/DESEQ_FILE.csv"
norm_counts_file <- "PATH/TO/norm_counts.csv"
columns <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8")
labels_columns <- c("CTL_1", "CTL_2", "CTL_3", "CTL_4",
               "EXP_5", "EXP_6", "EXP_7", "EXP_8")

# parameters
basemean_limit <- 100
up_log2FoldChange_limit <- 0.7
down_log2FoldChange_limit <- -0.7
padj_limit <- 0.05


info_data <- read.csv(deseq_file, header = T, row.names = 1)

info_data$log2FoldChange <- info_data$log2FoldChange * -1

norm_counts <- read.csv(norm_counts_file, header = T, row.names = 1)
total_reads <- colSums(norm_counts)
norm_counts <- norm_counts/total_reads*1000000
norm_counts <- log2(norm_counts)

upreg <- subset(info_data,
                baseMean > basemean_limit 
                & log2FoldChange > up_log2FoldChange_limit 
                & padj < padj_limit)
norm_counts$names <- rownames(norm_counts)
upreg$names <- rownames(upreg)
up_norm <- subset(norm_counts,
                  norm_counts$names %in% upreg$names,
                  select = columns)
up_norm <- as.matrix(up_norm)

downreg <- subset(info_data,
                  baseMean > basemean_limit 
                  & log2FoldChange > down_log2FoldChange_limit 
                  & padj < padj_limit)
downreg$names <- rownames(downreg)
down_norm <- subset(norm_counts,
                    norm_counts$names %in% downreg$names,
                    select = columns)
down_norm <- as.matrix(down_norm)

library(pheatmap)
pheatmap(up_norm, scale = "none",
         cellwidth = 20, cellheight = 15,
         cluster_rows = TRUE, cluster_cols = FALSE,
         labels_col = labels_columns)
pheatmap(down_norm, scale = "none",
         cellwidth = 20, cellheight = 10,
         cluster_rows = TRUE, cluster_cols = FALSE,
         labels_col = labels_columns)


