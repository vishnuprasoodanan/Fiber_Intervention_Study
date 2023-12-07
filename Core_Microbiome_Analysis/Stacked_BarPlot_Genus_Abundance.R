#Stacked bar plot

library(ggplot2)
library(tidyr)

df <- read.table(file="Mice.txt", sep = "\t", header = T)
data_order <- df[order(apply(df, 1, min)),]
# convert data from wide to long format using tidyr
df_long <- pivot_longer(df, cols = -Name, names_to = "Sample", values_to = "Abundance")
#colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "deepskyblue", "goldenrod1", "darkolivegreen", "salmon", "mediumpurple3", "paleturquoise4", "goldenrod4", "#FF00FF", "#F433FF", "#E600FF", "#CC00FF", "#B266FF","#9B30FF", "#8800FF", "#660066", "#993399", "#CC66CC", "#FF33CC", "#FF0099", "gray")
colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "gray")
# create stacked barplot using ggplot2
ggplot(df_long, aes(x = Sample, y = Abundance, fill = Name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "Abundance of Features in Each Sample", x = "Sample", y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

