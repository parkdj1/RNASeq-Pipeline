# ###########################
#
#  deseq2.R
#
#  Created on Wed May 6 2020
#  @author: parkdj1
#
#  run DeSeq2 on Salmon quant files
#
#  input: RNASeq quantification files
#  output files: normalized counts CSV file | paired reads CSV file
#  output plots: PCA plot | heatmap | MA plot
#
# ###########################

# set working directory (AKA where are your files going to be/output to)
setwd('../')

#change these to match data
o_folder = "PATH/TO/OUTPUT_FOLDER"
file_name_before_num = "PATH/TO/FILES"
file_name_after_num = "quant.sf"
column_data_file = "PATH/TO/column_data.txt"
num_samples = 10
CTL = "CTL_NAME"
EXP = "EXP_NAME"

list_samps = c(1:num_samples)
condition_values = c(rep(CTL, num_samples/2), rep(EXP, num_samples/2))


library(DESeq2)
library(ggplot2)
library("RColorBrewer")
library(pheatmap)

# process files and run through DeSeq2

# combine all modified quant files (Name/NumCounts cols only)
# files modified through python to have int counts and no 0's
my_data <-
  read.table(paste(file_name_before_num, "1", file_name_after_num, sep = ""),
             header = TRUE, sep = "\t")
colnames(my_data)[-1] <- "1"

for (num in list_samps[-1]) {
  file_name <- paste(file_name_before_num, num, file_name_after_num, sep = "")
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
my_data_filtered[, num_samples + 1] <- NULL

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

# PCA PLOT
plotPCA(rld, intgroup = c("Condition1", "Type"))
#   SAVE PLOT   #

rld_sampledist <- dist(t(assay(rld)))
rld_sampledistmatrix <- as.matrix(rld_sampledist)
rownames(rld_sampledistmatrix) <-
  paste(rld$Condition1, rld$Type, sep = "-")
colnames(rld_sampledistmatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# HEAT MAP
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
write.csv(counts, file = paste(o_folder,"/norm_counts.csv", sep = ""))
dds <- DESeq(dds)

# MA PLOT
plotMA(dds, ylim = c(-10,10))
#   SAVE PLOT   #

ddsPaired <- results(dds, c("Condition1", CTL, EXP))
write.csv(ddsPaired, paste(o_folder,"/",CTL,"_v_",EXP,".csv", sep = ""))


