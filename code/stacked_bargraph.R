setwd("Desktop/PUA_LAB/Genboree_output_AA_mice_2019")
library(ggplot2)

df <- read.csv("/Users/deborah/Desktop/PUA_LAB/Genboree_output_AA_mice_2019/rpm_prop.csv")
head(df)

ggplot(df, aes(fill=YRNA, x=count, y=type, label=count)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x="Proportion of RPM for RNY1 and RNY3", y="")
ggsave("bargraph_unfilled_yrna_prop_rpm.png",width=5,height=3)
ggplot(df, aes(fill=YRNA, x=count, y=type)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x="Proportion of RPM for RNY1 and RNY3", y="")
ggsave("bargraph_filled_yrna_prop_rpm.png",width=5,height=3)
