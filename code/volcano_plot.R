# ###########################
#
#  volcano_plot.R
#
#  Created on Wed May 6 2020
#  @author: parkdj1
#
#  create a volcano plot from CSV data
#
#  input: paired reads CSV file
#  output: volcano plot
#
# ###########################

# set working directory (AKA where are your files going to be/output to)
setwd('../')

# set path to CSV file containing data
my_csv <- "PATH/TO/PAIRED_READS.csv"

# set Volcano graph parameters
graph_title <- "TITLE_OF_GRAPH"
x_data <- "X_COLUMN_NAME"
y_data <- "Y_COLUMN_NAME"
pCutoff <- 0.01
FCcutoff <- 1
plotname <- "OUTPUT_FILE_NAME.png"

data <- read.csv(my_csv, header = T, row.names = 1)

data$log2FoldChange <- data$log2FoldChange * -1
data <- subset(data, baseMean > 100)

library(EnhancedVolcano)

EnhancedVolcano(data,
                lab = rownames(data),
                x = x_data,
                y = y_data,
                xlim = c(-8, 8),
                ylim = c(-1,22),
                title = graph_title,
                pCutoff = pCutoff,
                FCcutoff = FCcutoff,
                pointSize = 3.0,
                labSize = 0,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)


