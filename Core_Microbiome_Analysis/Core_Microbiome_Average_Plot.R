#Plot Avg. relative abundance of Core microbiomes in Citrobacter FMT samples and initial human samples.
#There are three input files; (1) Relative abundance: Core microbiome of Citrobacter FMT samples, (2) Relative abundance: Core Microbiome of Human samples in visit-1
#(3) Relative abundance: Core microbiome of Human samples in visit-3
#The input file format: rows contain the genera and column contain sample names. The last row contain the Group information for each sample
library(ggplot2)
library(dplyr)
library(plyr)
set.seed(123)

df1 <- read.table("Citro_Core_80.txt", sep="\t", header=T)
df2 <- read.table("Visit_1_Core_80.txt", sep="\t", header=T)
df3 <- read.table("Visit_3_Core_80.txt", sep="\t", header=T)
# Merge dataframes based on common row names (intersection)
merged_df1 <-merge(df1, df2, by = "Genera", all = FALSE)
merged_df2 <-merge(merged_df1, df3, by = "Genera", all = FALSE)

names(merged_df2) <- gsub("^X", "", names(merged_df2))
names(merged_df2) <- gsub("\\.", "-", names(merged_df2))


### plot stacked bar-plots using the average relative abundance

data <- as.data.frame(t(merged_df2))
colnames(data) <- data[1,]
data <- data[-1, ]
class <- data$Status
data <- data %>% mutate_if(~ !identical(., last(.)) && is.character(.), as.numeric)
data$Status <- class

group_averages <- data %>% group_by(Status) %>% dplyr::summarize(across(starts_with("d__Bacteria"), mean))
group_averages_cs <- colSums(group_averages[, -1])

# Find the names of the columns with the 20 highest column sums
group_averages_top_columns <- names(group_averages_cs)[order(group_averages_cs, decreasing = TRUE)[1:5]]

# Select the columns with the 20 highest column sums
group_averages_top_30 <- group_averages[, c("Status", group_averages_top_columns)]
group_averages_top_31 <- group_averages_top_30 %>% mutate(Others = 1-(rowSums(select(., -1))))

write.table(group_averages_top_31, file = "Core25_Average_Values.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
data2 <- read.table("Core25_Average_Values_edited.txt", sep="\t", header=T)
df_long <- pivot_longer(data2, cols = -Status, names_to = "Sample", values_to = "Abundance")

colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "deepskyblue","gray")
#colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "deepskyblue", "goldenrod1", "darkolivegreen", "salmon", "mediumpurple3", "paleturquoise4", "goldenrod4", "#FF00FF", "#F433FF", "#E600FF", "#CC00FF", "#B266FF","#9B30FF", "#8800FF", "#660066", "#993399", "#CC66CC", "#FF33CC", "#FF0099", "tomato", "limegreen", "dodgerblue", "purple", "orchid", "saddlebrown", "turquoise", "gray")

# create stacked barplot using ggplot2
ggplot(df_long, aes(x = Status, y = Abundance, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "Core Mirobiome (80%)", x = "Status", y = "Average Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
