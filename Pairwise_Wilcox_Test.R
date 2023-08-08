##---------------------------------------------------------
# wilcoxon test for Human samples
# location: /Users/vishnu/WORK/FIBER_INTERVENTION_16S/CORE_MICROBIOME_STACKED_BARPLOTS/WICOX_TEST

# initialize an empty data frame
results_df <- data.frame()
V1 <- read.table(file = "human_v1.txt", sep = "\t", header = T, row.names = 1)
V3 <- read.table(file = "human_v3.txt", sep = "\t", header = T, row.names = 1)
for (col_name in colnames(V1)) {
  name <- col_name
  pvalue <- wilcox.test(V1[[col_name]], V3[[col_name]], alternative = "two.sided", paired = TRUE)$p.value
  # create a new row of data using rbind
  new_row <- data.frame(name, pvalue)
  
  # add the new row to the results data frame using rbind
  results_df <- rbind(results_df, new_row)
}
write.table (results_df, file="Wilcox_test.txt", sep = "\t", row.names = FALSE, quote = FALSE)
