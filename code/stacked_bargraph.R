# ###########################
#
#  stacked_bargraph.R
#
#  Created on Wed May 6 2020
#  @author: parkdj1
#
#  create stacked bar graphs (filled/unfilled) from CSV data
#
#  input: paired reads CSV file
#  output: filled stacked bar graph | unfilled stacked bar graph
#
# ###########################

# set working directory (AKA where are your files going to be/output to)
setwd('../')

# set path to CSV file containing data
my_csv <- "PATH/TO/PAIRED_READS.csv"

# set plot data
unfilled_plot_name <- paste("PLOTNAME", "unfilled.png", sep = "_")
unfilled_plot_x_label <- "Proportion of RPM for RNY1 and RNY3"
unfilled_fill_col <- "YRNA"
unfilled_x_col <- "count"
unfilled_y_col <- "type"

filled_plot_name <- paste("PLOTNAME", "filled.png", sep = "_")
filled_plot_x_label <- "Proportion of RPM for RNY1 and RNY3"
filled_fill_col <- "YRNA"
filled_x_col <- "count"
filled_y_col <- "type"

library(ggplot2)

df <- read.csv(my_csv)
head(df)

# unfilled bar graph
ggplot(df, aes(fill = unfilled_fill_col,
               x = unfilled_x_col,
               y = unfilled_y_col,
               label = unfilled_x_col)) + 
  geom_bar(position = "stack", stat = "identity") +
  labs(x = unfilled_plot_x_label, y = "")
ggsave(unfilled_plot_name, width = 5, height = 3)

# filled bar graph
ggplot(df, aes(fill = filled_fill_col,
               x = filled_x_col,
               y = filled_y_col)) + 
  geom_bar(position = "fill", stat = "identity") +
  labs(x = filled_plot_x_label, y = "")
ggsave(filled_plot_name, width = 5, height = 3)

