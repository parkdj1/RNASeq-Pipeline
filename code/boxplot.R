# ###########################
#
#  boxplot.R
#
#  Created on Wed May 6 2020
#  @author: parkdj1
#
#  create box plots from CSV data
#
#  input: CSV file
#  output: box plot
#
# ###########################

# set working directory (AKA where are your files going to be/output to)
setwd('../')

# set path to CSV file containing data
my_csv <- "PATH/TO/DATAFILE.csv"

library(ggplot2)
library(outliers)

df <- read.csv(my_csv)
head(df)
grubbs.test(df$count)
df_no_outlier <- head(df,-1)
df1 <- subset(df,YRNA=='RNY1', select=c(type,count))
df3 <- subset(df,YRNA=='RNY3', select=c(type,count))

p <- ggplot(df, aes(x=type, y=count,fill=YRNA))+
  geom_boxplot()+
  coord_flip()+
  labs(y="Y RNA Reads Per Million", x="")+
  facet_wrap(~YRNA, scales="free_x")
p
p1 <- ggplot(df1, aes(x=type, y=count))+
  geom_boxplot(color="black",fill="#F8766D")+
  coord_flip()+
  labs(y="Y RNA Reads Per Million", x="")
p1
p3 <- ggplot(df3, aes(x=type, y=count))+
  geom_boxplot(color="black",fill="#00BFC4")+
  coord_flip()+
  labs(y="Y RNA Reads Per Million", x="")
p3

p_no_outlier <- ggplot(df_no_outlier, aes(x=type, y=count,fill=YRNA))+
  geom_boxplot()+
  coord_flip()+
  labs(y="Y RNA Reads Per Million", x="")+
  facet_wrap(~YRNA, scales="free_x")
p_no_outlier
ggsave("boxplot_rpm_prop-excl-outlier.png",width=5,height=3)

p
ggsave("boxplot_rpm_prop.png",width=5,height=3)
p1
ggsave("YRNA1_boxplot_rpm_prop.png",width=5,height=3)
p3
ggsave("YRNA3_boxplot_rpm_prop.png",width=5,height=3)


