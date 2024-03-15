#This script necessitates three input files. The 'Human_Genus_count_table.txt' file comprises genus (or ASV) IDs listed in rows and the total count of assigned genera/ASVs in each sample denoted in columns. 'Taxonomy.tsv' includes the genus/ASV IDs in the first column and their corresponding taxonomy in the second column. Lastly, 'Humans_Metadata.txt' contains sample names in the first column and their respective metadata (e.g., individual ID, HD or HF, etc.) listed in the remaining columns.

# create PCoA using clr-transformed data
counts  <- read.csv("Human_Genus_count_table.txt", sep = "\t", row.names = 1)  # Abundance table (e.g. ASV data; to assay data)
tax     <- as.matrix(read.csv("taxonomy.tsv", sep = "\t", row.names = 1))     # Taxonomy table (to rowData)
samples <- as.matrix(read.table("Humans_Metadata.txt", sep="\t", header = TRUE, row.names = 1))  # collate data (to colData)
#counts <- as.data.frame(counts[, -1])
names(counts) <- gsub("^X", "", names(counts))
names(counts) <- gsub("\\.", "-", names(counts))
counts[, -1] <- apply(counts[, -1], 2, function(x) as.numeric(x))
#counts <- apply(counts, 2, function(x) as.numeric(x))
counts <- as.matrix(counts)  
se <- SummarizedExperiment(assays = list(counts = counts),
                           colData = samples,
                           rowData = tax)
tse <- as(se, "TreeSummarizedExperiment")
tse <- transformAssay(se, method = "clr", pseudocount = 1)
clr_assay <- assays(tse)$clr
clr_assay <- t(clr_assay)

# Calculates Euclidean distances between samples. Because taxa is in columns,
# it is used to compare different samples.
euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")
Euclidean_distances <- as.matrix(euclidean_dist)
Euclidean_dist_column <- melt(Euclidean_distances)
write.table (Euclidean_dist_column, file = "Euclidean_distances.txt", sep = "\t")

eucl_pcoa <- pcoa(Euclidean_distances)
eucl_pcoa$values[1:2,]
mds.var.per = round(eucl_pcoa$values$Eigenvalues/sum(eucl_pcoa$values$Eigenvalues)*100, 1)
pc <- c(1,2)
pcoa_df <- eucl_pcoa$vectors[,1:2]
jpeg("PCoA.jpg", height = 10, width = 10, units = 'in', res = 600)
plot(eucl_pcoa$vectors[,1:2], bg= c("skyblue","palegreen", "orangered")[as.factor(as.data.frame(samples)$Status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
ordiellipse(eucl_pcoa$vectors[,1:2], as.data.frame(samples)$Status, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("skyblue","palegreen", "orangered"))
#legend("bottomleft", legend = c("HD", "HF"), col = c("skyblue","palegreen", "orangered"), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")

for (i in 1:(length(as.data.frame(samples)$Status3)-1)) {
  if (!is.na(as.data.frame(samples)$Status3[i])) {
    for (j in (i+1):length(as.data.frame(samples)$Status3)) {
      if (!is.na(as.data.frame(samples)$Status3[j]) && as.data.frame(samples)$Status3[i] == as.data.frame(samples)$Status3[j]) {
        x <- c(pcoa_df[i, "Axis.1"], pcoa_df[j, "Axis.1"])
        y <- c(pcoa_df[i, "Axis.2"], pcoa_df[j, "Axis.2"])
        lines(x, y, col = "black", lty = "solid")
      }
    }
  }
}
dev.off ()
adonis2(Euclidean_distances ~ as.data.frame(samples)$Status)
