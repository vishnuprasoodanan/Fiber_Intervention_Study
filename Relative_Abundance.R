df <- read.table(file="genera_abundance.txt", sep = "\t", header = T, row.names = 1)
names(df) <- gsub("^X", "", names(df))
names(df) <- gsub("\\.", "-", names(df))
# calculate the relative abundance
df <- as.data.frame(t(df))
df_relative <- t(apply(df, 1, function(x) x/sum(x)))

# transpose the result back to its original orientation
df_relative <- as.data.frame(t(df_relative))
write.table(df_relative, file = "Genus_relative_ab.txt", sep = '\t', quote = FALSE)
