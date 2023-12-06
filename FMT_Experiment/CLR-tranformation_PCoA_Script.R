# create PCoA using clr-transformed data
counts  <- read.csv("selected-feature-table.txt", sep = "\t", row.names = 1)  # Abundance table (e.g. ASV data; to assay data)
tax     <- as.matrix(read.csv("selected_taxonomy.txt", sep = "\t", row.names = 1))     # Taxonomy table (to rowData)
samples <- as.matrix(read.table("Metadata_selected.tsv", sep="\t", header = TRUE, row.names = 1))  # collate data (to colData)
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

# Does principal coordinate analysis
# euclidean_pcoa <- ecodist::pco(euclidean_dist)
# 
#  # Creates a data frame from principal coordinates
#  euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1],
#                                  pcoa2 = euclidean_pcoa$vectors[,2])
# 
#  euclidean_patient_status_pcoa_df <- cbind(euclidean_pcoa_df,
#                                            patient_status = colData(tse)$Status)
# 
#  # Creates a plot
#  euclidean_patient_status_plot <- ggplot(data = euclidean_patient_status_pcoa_df,
#                                          aes(x=pcoa1, y=pcoa2,
#                                              color = patient_status)) +
#    geom_point() +
#    labs(x = "PC1",
#         y = "PC2",
#         title = "PCoA with Aitchison distances") +
#    theme(title = element_text(size = 12)) # makes titles smaller
# 
#  euclidean_patient_status_plot

#UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
UniFrac_distances <- as.matrix(euclidean_dist)
UniFrac_dist_column <- melt(UniFrac_distances)
write.table (UniFrac_dist_column, file = "UniFrac_distances.txt", sep = "\t")

Uni_pcoa <- pcoa(UniFrac_distances)
Uni_pcoa$values[1:2,]
mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
pc <- c(1,2)

jpeg("UniFrac_PCoA.jpg", height = 10, width = 10, units = 'in', res = 600)
plot(Uni_pcoa$vectors[,1:2], bg= c("darkolivegreen4", "red", "darkgreen", "salmon4")[as.factor(as.data.frame(samples)$Status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
ordiellipse(Uni_pcoa$vectors[,1:2], as.data.frame(samples)$Status, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("darkolivegreen4", "red", "darkgreen", "salmon4"))
#ordispider(Uni_pcoa$vectors[,1:2], as.data.frame(samples)$Status, lty=3, spider ="centroid", lwd=1, col="black")
#legend("bottomleft", legend = c("A_BL_Chow", "B_BL_WSD", "C_HF_Chow", "D_HF_WSD"), col = c("darkolivegreen4", "red", "darkgreen", "salmon4"), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
adonis2(UniFrac_distances ~ as.data.frame(samples)$Status)
