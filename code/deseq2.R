# ############
# Deborah Park 05.06.2020
# ###########

library(DESeq2)
library(ggplot2)

#change these to match data
setwd("Desktop/PUA_LAB")
folder = "Serezani"
file_before_num = "Serezani/salmon_quant/"
file_after_num = "quant_collapsed.sf"
column_data_file = "Serezani/column_data.txt"
num_samples = 10
list_samps = c(1:10)
CTL = "NOT_INFECTED"
EXP = "mkSA"
condition_values = c(rep(CTL, 5), rep(EXP, 5))

# process files and run through DeSeq2

# combine all modified quant files (Name/NumCounts cols only)
# files modified through python to have int counts and no 0's
my_data <-
  read.table(paste(file_before_num, "1", file_after_num, sep = ""), header = TRUE, sep = "\t")
colnames(my_data)[-1] <- "1"

for (num in list_samps[-1]) {
  file_name <- paste(file_before_num, num, file_after_num, sep = "")
  new_col <- read.table(file_name, header = TRUE, sep = "\t")
  colnames(new_col)[-1] <- num
  my_data <- merge(my_data, new_col, by = "Name", all = TRUE)
}

# write all null to 1
my_data[is.null(my_data)] <- 1
# change column to row names
rownames(my_data) <- my_data[, 1]
my_data[, 1] <- NULL
# filter low-expressed genes
my_data$averages <- rowMeans(my_data)
my_data_filtered <- subset(my_data, averages > 10)
my_data_filtered[, num_samples+1] <- NULL

# other files necessary for DeSeq
column_data <-
  read.table(column_data_file, header = T, row.names = 1)
condition <-
  DataFrame(factor(condition_values))

# run DESeq
dds <-
  DESeqDataSetFromMatrix(countData = my_data_filtered,
                         colData = column_data,
                         design = ~ Condition1)
rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup = c("Condition1", "Type"))
#   SAVE PLOT   #

rld_sampledist <- dist(t(assay(rld)))
library("RColorBrewer")
library(pheatmap)
rld_sampledistmatrix <- as.matrix(rld_sampledist)
rownames(rld_sampledistmatrix) <-
  paste(rld$Condition1, rld$Type, sep = "-")
colnames(rld_sampledistmatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(
  rld_sampledistmatrix,
  clustering_distance_rows = rld_sampledist,
  clustering_distance_cols = rld_sampledist,
  col = colors
)
#   SAVE PLOT   #

#Get back normalized count values from DESeq2
dds <- estimateSizeFactors(dds)
counts <- counts(dds, normalized = TRUE)
write.csv(counts, file = paste(folder,"/norm_counts.csv",sep=""))
dds <- DESeq(dds)
plotMA(dds, ylim = c(-10,10))
#   SAVE PLOT   #
ddsPaired <- results(dds, c("Condition1", CTL, EXP))
write.csv(ddsPaired, paste(folder,"/",CTL,"_v_",EXP,".csv",sep=""))

