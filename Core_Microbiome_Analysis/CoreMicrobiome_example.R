#This script will evaluate core microbiome with a criteria that the genus should be present in 80% of the samples in the group
library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

df <- read.table(file="Genus_count_table.txt", sep = "\t", header = T, row.names = 1)
names(df) <- gsub("^X", "", names(df))
names(df) <- gsub("\\.", "-", names(df))
# calculate the relative abundance
df <- as.data.frame(t(df))
df_relative <- t(apply(df, 1, function(x) x/sum(x)))

#-----------------Generating core microbiome for FMT mice samples
#read the 'Genus_count_table.txt' file. Nemas of genera in rows and samplesnames in columns 
#data <- read.table(file = "Genus_count_table.txt", sep = "\t", header = T, row.names = 1)
data_t <- as.data.frame(df_relative)
rownames(data_t) <- gsub("\\.", "_", rownames(data_t))
# Subset rows with rownames starting with 'HF-Chow', 'HD-WSD', 'HF-Chow', 'HF-WSD' and exclude rownames ending with '-BL'

FMT_mice_data <- data_t %>%
  rownames_to_column(var = "rowname") %>%
  filter(grepl("^(HF-Chow|HD-Chow|HD-WSD|HF-WSD)", rowname) & !grepl("-BL$", rowname))

# Remove row names and convert the first column to row names
rownames(FMT_mice_data) <- FMT_mice_data[, 1]
FMT_mice_data <- FMT_mice_data[, -1]

FMT_mice_data <- FMT_mice_data[, colSums(FMT_mice_data) != 0.0]

# Calculate 80% of the number of rows
num_rows <- nrow(FMT_mice_data) 
eighty_cutoff <- as.numeric(round(0.8 * num_rows, 0.0))  

# Identify columns with non-zero values in 'eighty_cutoff' or more rows
selected_columns <- which(colSums(FMT_mice_data != 0.0) >= eighty_cutoff)

# Subset the dataframe to keep only the selected columns
FMT_80_core <- FMT_mice_data[, selected_columns]

#-------------Generating core microbiome for Human samples with habitual diet
human_hd <- data_t %>%
  rownames_to_column(var = "rowname") %>%
  filter(grepl("^HD_0", rowname))

# Remove row names and convert the first column to row names
rownames(human_hd) <- human_hd[, 1]
human_hd <- human_hd[, -1]

human_hd <- human_hd[, colSums(human_hd) != 0.0]

# Calculate 80% of the number of rows
num_rows <- nrow(human_hd) 
eighty_cutoff <- as.numeric(round(0.8 * num_rows, 0.0))  

# Identify columns with non-zero values in 'eighty_cutoff' or more rows
selected_columns <- which(colSums(human_hd != 0.0) >= eighty_cutoff)

# Subset the dataframe to keep only the selected columns
human_hd_80_core <- human_hd[, selected_columns]


#Generating core microbiome for Human samples with high fiber diet
human_hf <- data_t %>%
  rownames_to_column(var = "rowname") %>%
  filter(grepl("^HF_0", rowname))

# Remove row names and convert the first column to row names
rownames(human_hf) <- human_hf[, 1]
human_hf <- human_hf[, -1]

human_hf <- human_hf[, colSums(human_hf) != 0.0]

# Calculate 80% of the number of rows
num_rows <- nrow(human_hf) 
eighty_cutoff <- as.numeric(round(0.8 * num_rows, 0.0))  

# Identify columns with non-zero values in 'eighty_cutoff' or more rows
selected_columns <- which(colSums(human_hf != 0.0) >= eighty_cutoff)

# Subset the dataframe to keep only the selected columns
human_hf_80_core <- human_hf[, selected_columns]

# Get column names of each dataframe
col_names_df1 <- colnames(FMT_80_core)
col_names_df2 <- colnames(human_hd_80_core)
col_names_df3 <- colnames(human_hf_80_core)

# Find common genera
common_columns <- intersect(intersect(col_names_df1, col_names_df2), col_names_df3)

# Display common genera
print(common_columns)

common_FMT_mice_80_core <- FMT_80_core[, common_columns]
common_human_hd_80_core <- human_hd_80_core[, common_columns]
common_human_hf_80_core <- human_hf_80_core[, common_columns]

# Add dataframe names as the first column
common_FMT_mice_80_core$Status <- "FMT_mice"
common_human_hd_80_core$Status <- "human_hd"
common_human_hf_80_core$Status <- "human_hf"

# Merge the dataframes by rows and add dataframe names as the first column
merged_df <- bind_rows(common_FMT_mice_80_core, common_human_hd_80_core, common_human_hf_80_core)

# Reorder columns (optional)
merged_df <- merged_df[, c(ncol(merged_df), 1:(ncol(merged_df) - 1))]

sum_columns_2_to_5 <- rowSums(merged_df[, 2:5])

# Calculate 1 - sum of columns 2, 3, 4, and 5
merged_df$Others <- 1 - sum_columns_2_to_5
merged_df_v1 <- merged_df

#----------------------------------------------------------------------
#sample wise abundance plot
merged_df <- merged_df[, !colnames(merged_df) %in% 'Status', drop = FALSE]
merged_df <- as.data.frame(t(merged_df))

merged_df$Name <- row.names(merged_df)

# Reset row names to NULL to remove them from the dataframe
row.names(merged_df) <- NULL

# Reorder columns (optional if you want the newly added column as the first column)
merged_df <- merged_df[, c(ncol(merged_df), 1:(ncol(merged_df) - 1))]

data_order <- merged_df[order(apply(merged_df, 1, min)),]
# convert data from wide to long format using tidyr

df_long <- pivot_longer(merged_df, cols = -Name, names_to = "Sample", values_to = "Abundance")
colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "gray")
# create stacked barplot using ggplot2
ggplot(df_long, aes(x = Sample, y = Abundance, fill = Name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "Abundance of Features in Each Sample", x = "Sample", y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

###--------------------------------------------------------------------------
### plot stacked bar-plots using the average relative abundance
library(ggplot2)
library(dplyr)
set.seed(123)
data1 <- merged_df_v1
# Save column names into a vector
column_names <- colnames(data1)
column_names <- column_names[-1]

group_averages<- data1 %>% group_by(Status) %>% dplyr::summarize(across(column_names, mean))
write.table(group_averages, file = "group_average_rel_abundance.txt", sep = '\t', quote = FALSE, row.names = FALSE)

df2_long <- pivot_longer(group_averages, cols = -Status, names_to = "Sample", values_to = "Abundance")
colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "gray")
# create stacked barplot using ggplot2
ggplot(df2_long, aes(x = Status, y = Abundance, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "Abundance of Features in Each Sample", x = "Status", y = "Average Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
