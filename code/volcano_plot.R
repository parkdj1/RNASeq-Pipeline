setwd('../')

data <- read.csv("neil/deseq2_collapsed/HFD_ff_v_HFD_ko.csv",header=T,row.names=1)

data$log2FoldChange <- data$log2FoldChange * -1
data <- subset(data, baseMean > 100)

library(EnhancedVolcano)

EnhancedVolcano(data,
                lab = rownames(data),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-8, 8),
                ylim = c(-1,22),
                title = 'all samples, by gene id, cutoff 10, basemean >100',
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 0,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)
