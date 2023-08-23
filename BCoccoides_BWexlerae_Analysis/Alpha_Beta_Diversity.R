library(phyloseq)
library (ape)
library (ggplot2)
library("GUniFrac")
#----------------------------------- Alpha diversity analysis on complete data (Total-Data)
#-- Make a phyloseq object after filtering of OTUs
# Read in metadata
metadata <- read.table("Metadata.tsv", sep="\t", header=T)
# Subset "SampleID" based on "treatment" condition
#selected_sample_ids <- metadata$SampleID[metadata$Treatment %in% c("BT_HK", "BT_live", "Control", "BW_HK", "BW_live")]
selected_sample_ids <- metadata$SampleID[metadata$Treatment %in% c("BT_HK", "BT_live", "Control")]
# Exclude the selected SampleID entries from the dataframe
filtered_metadata <- subset(metadata, !(SampleID %in% selected_sample_ids))
# Set the first column as rownames
rownames(filtered_metadata ) <- filtered_metadata$SampleID
# Remove the first column
filtered_metadata <- filtered_metadata[, -1]
filtered_metadata$Status <- paste(filtered_metadata$Timepoint, filtered_metadata$Treatment, sep = "_")

# Read in OTU table
otu_table_in <- read.csv("feature-table.tsv", sep = "\t", row.names = 1)
names(otu_table_in) <- gsub("^X", "", names(otu_table_in))
names(otu_table_in) <- gsub("\\.", "-", names(otu_table_in))
# Exclude columns with names present in selected_sample_ids
otu_table_filtered_data <- otu_table_in[, !(colnames(otu_table_in) %in% selected_sample_ids)]
# Exclude rows with row sum equal to zero
otu_table_filtered_data_final <- otu_table_filtered_data[rowSums(otu_table_filtered_data) != 0, ]

#OTU_table <- read.table(file = "selected_feature-table.txt", sep = "\t", header = T, row.names = 1) #ASV-IDs in rows and Samples in columns
#OTU_table <- as.data.frame(t(otu_table_filtered_data_final)) #Samples in rows and ASV-IDs in columns
# <- Rarefy(OTU_table, depth = min(rowSums(OTU_table)))$otu.tab.rff # Samples in Row and OTUs in column
#OTU_table_rar <- data.frame(otu.tab.rff)

#write.table(OTU_table_rar, file = "Selected_FeatureTable_rarefied.txt", quote = FALSE, sep = '\t')
#OTU_table_rar_final <- as.data.frame(t(OTU_table_rar))

total_asv <- colSums(otu_table_filtered_data_final)
range (total_asv)
sd (total_asv)
mean (total_asv)

# Read in taxonomy
# Separated by kingdom, phylum, class, order, family, genus, species
taxonomy <- read.csv("taxonomy.tsv", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)

# Read in tree
phy_tree <- read_tree("tree.nwk")
# Import all as phyloseq objects
OTU <- otu_table(otu_table_filtered_data_final, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(filtered_metadata)
# Sanity checks for consistent OTU names
head(taxa_names(TAX))
head(taxa_names(OTU))
head(taxa_names(phy_tree))
# Same sample names
sample_names(OTU) <- gsub('[.]', '-', sample_names(OTU))
sample_names(OTU) <- gsub('^[X]', '', sample_names(OTU)) 
head(sample_names(OTU))
head(sample_names(META))

# Finally merge to create Phyloseq object!
ps <- phyloseq(OTU, TAX, META, phy_tree)

covar <- as.vector(colnames(filtered_metadata))
# Elements to remove
elements_to_remove <- c("Mouse_number", "Sample_type")

# Remove specific elements
covar_upd <- covar[!(covar %in% elements_to_remove)]
for (i in 1:length(covar_upd)) {
  element <- covar_upd[i]
  if (element == "Timepoint") {
    Colors <- c("darkolivegreen4", "red")
    Colors1 <- c("darkolivegreen","red4")
  } else if (element == "Treatment"){
    Colors <- c("mediumpurple3", "snow4", "aquamarine", "orange")
    Colors1 <- c("violet", "gray30", "deepskyblue", "orangered")
  }else if (element == "Status"){
    Colors <- c("wheat", "goldenrod1", "mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine", "salmon4")
    Colors1 <- c("wheat4", "goldenrod4", "mediumpurple4", "grey27", "blue", "paleturquoise4", "aquamarine4",  "salmon")
  }else if (element == "Sex"){
    Colors <- c("lightblue", "gray")
    Colors1 <- c("midnightblue", "darkgray")
  }
  #f_col1 <- paste(element, "Colors1", sep = "_")
  #f_col <- paste(element, "Colors", sep = "_")
  p <- plot_richness(ps, element, measures = c("Observed", "Chao1", "Shannon"), color = NULL, shape = NULL)
  p <- p + geom_boxplot(data = p$data, aes_string(x= element, y = p$data$value, fill = element)) + 
    labs(x="",y="Alpha Diversity Measure") + 
    theme_bw() +
    scale_color_manual(values = Colors1)+
    scale_fill_manual(values = Colors)+
    geom_point(aes_string(colour = element), position=position_jitterdodge(jitter.width = 0.5), size = 3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 24))
  file_name1 <- paste(element, "phyloseq_analysis-richness_estimates.pdf", sep = "_")
  file_name2 <- paste(element, "phyloseq_analysis-richness_estimates.jpg", sep = "_")
  ggsave(file_name1, p, width = 10, height = 10, limitsize = FALSE)
  jpeg(file_name2, height = 10, width = 15, units = 'in', res = 600)
  p
  dev.off ()
  #---------------------------PCoA Analysis
  
  Bushman2  = transform_sample_counts(ps, function(x) x / sum(x) )
  
  #---------------------- weighted UniFrac PCoA
  file_name3 <- paste(element, "W-UniFrac_distances.txt", sep = "_")
  UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
  UniFrac_distances <- as.matrix(UniFrac_distances)
  UniFrac_dist_column <- melt(UniFrac_distances)
  write.table (UniFrac_dist_column, file = file_name3, sep = "\t")
  
  Uni_pcoa <- pcoa(UniFrac_distances)
  Uni_pcoa$values[1:2,]
  mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
  pc <- c(1,2)
  file_name4 <- paste(element, "W-Unifrac_PCOA.jpg", sep = "_")
  jpeg(file_name4, height = 10, width = 10, units = 'in', res = 600)
  plot(Uni_pcoa$vectors[,1:2], bg= Colors1[as.factor(Bushman2@sam_data[[element]])], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
  text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
  ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = Colors1)
  #ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], lty=3, spider ="centroid", lwd=1, col="black")
  legend("bottomleft", legend = sort(unique(as.factor(Bushman2@sam_data[[element]]))), col = Colors1,lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
  abline(h=0, v=0, col = "gray60")
  dev.off ()
  #adonis2(UniFrac_distances ~ Bushman2@sam_data[[element]])
  
  #--------------------- unweighted UniFrac PCoA
  
  UWUniFrac_distances <- UniFrac(Bushman2, weighted=FALSE)
  UWUniFrac_distances <- as.matrix(UWUniFrac_distances)
  UniFrac_dist_column <- melt(UWUniFrac_distances)
  file_name5 <- paste(element, "unweighted_UWUniFrac_distances.txt", sep = "_")
  write.table (UniFrac_dist_column, file = file_name5, sep = "\t")
  Uni_pcoa <- pcoa(UWUniFrac_distances)
  Uni_pcoa$values[1:2,]
  mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
  pc <- c(1,2)
  
  file_name6 <- paste(element, "unweighted_UWUniFrac_distances.jpg", sep = "_")
  jpeg(file_name6, height = 10, width = 10, units = 'in', res = 600)
  plot(Uni_pcoa$vectors[,1:2], bg= Colors1[as.factor(Bushman2@sam_data[[element]])], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
  text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
  ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = Colors1)
  #ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], lty=3, spider ="centroid", lwd=1, col="black")
  legend("bottomleft", legend = sort(unique(as.factor(Bushman2@sam_data[[element]]))), col = Colors1,lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
  abline(h=0, v=0, col = "gray60")
  dev.off ()
  #adonis2(UWUniFrac_distances ~ Bushman2@sam_data$diet)
  
  #-------------------------- Bray-Curtis PCoA
  BC_distances <- distance(Bushman2, method="bray")
  BC_distances <- as.matrix(BC_distances)
  UniFrac_dist_column <- melt(BC_distances)
  
  file_name7 <- paste(element, "BC_distances.txt", sep = "_")
  write.table (UniFrac_dist_column, file = file_name7, sep = "\t")
  Uni_pcoa <- pcoa(BC_distances)
  Uni_pcoa$values[1:2,]
  mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
  pc <- c(1,2)
  
  file_name8 <- paste(element, "BC_distances.jpg", sep = "_")
  jpeg(file_name8, height = 10, width = 10, units = 'in', res = 600)
  plot(Uni_pcoa$vectors[,1:2], bg= Colors1[as.factor(Bushman2@sam_data[[element]])], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
  text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
  ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = Colors1)
  #ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data[[element]], lty=3, spider ="centroid", lwd=1, col="black")
  legend("topleft", legend = sort(unique(as.factor(Bushman2@sam_data[[element]]))), col = Colors1,lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
  abline(h=0, v=0, col = "gray60")
  dev.off ()
  #adonis2(BC_distances ~ Bushman2@sam_data$diet)
}
