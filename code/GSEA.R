# ###########################
#
#  GSEA.R
#
#  Created on Wed May 6 2020
#  @author: parkdj1
#
#  prepare deseq files for GSEA
#
#  input: normalized counts CSV file
#  output: modified normalized counts TXT file
#
# ###########################

# set working directory (AKA where are your files going to be/output to)
setwd('../')

# set path to CSV files containing data
norm_counts_file <- "PATH/TO/norm_counts.csv"
transcript_biotype <- "protein_coding"
columns <- c("X9", "X10", "X11", "X12", "X13", "X14", "X15", "X16")
col_names <- c("CTL_1", "CTL_2", "CTL_3", "CTL_4",
                   "EXP_5", "EXP_6", "EXP_7", "EXP_8")
avg_reads_limit <- 100
table_name <- "TABLE_NAME.txt"


norm_counts <- read.csv(norm_counts_file, header = T, row.names = 1)
norm_counts <- subset(norm_counts,
                      transcript_biotype == transcript_biotype)
norm_counts_num <- subset(norm_counts,
                          transcript_biotype == transcript_biotype,
                          select = columns)

norm_counts_num$avg_reads <- rowSums(norm_counts_num)/8
norm_counts_num <- subset(norm_counts_num,
                          avg_reads > avg_reads_limit,
                          select = columns)

norm_counts$names <- rownames(norm_counts)
norm_counts_num$names <- rownames(norm_counts_num)

norm_counts <- subset(norm_counts,
                      norm_counts$names %in% norm_counts_num$names,
                      select = c(gene_symbol, description, columns))

rownames(norm_counts) <- NULL
colnames(norm_counts) <- c("NAME","DESCRIPTION", col_names)

write.table(norm_counts, table_name, row.names = FALSE, sep = "\t")

