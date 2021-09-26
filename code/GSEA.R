setwd(file.path("Desktop", "PUA_LAB","deseq_4_excl_10_cutoff_10"))

norm_counts <- read.csv("norm_counts_mod.csv",header=T,row.names=1)
norm_counts <- subset(norm_counts, transcript_biotype == "protein_coding")
norm_counts_num <- subset(norm_counts, transcript_biotype == "protein_coding", select = c(X9,X11,X12,X13,X14,X15,X16))

norm_counts_num$avg_reads <- rowSums(norm_counts_num)/8
norm_counts_num <- subset(norm_counts_num, avg_reads > 100, select = c(X9,X11,X12,X13,X14,X15,X16))

norm_counts$names <- rownames(norm_counts)
norm_counts_num$names <- rownames(norm_counts_num)

norm_counts <- subset(norm_counts, norm_counts$names %in% norm_counts_num$names, select = c(gene_symbol,description,X9,X11,X12,X13,X14,X15,X16))

rownames(norm_counts) <- NULL
colnames(norm_counts) <- c("NAME","DESCRIPTION","ctr_9","ctr_11","ctr_12","cnot6_13","cnot6_14","cnot6_15","cnot6_16")

write.table(norm_counts, "GSEA_norm_counts.txt", row.names = FALSE, sep = "\t")
