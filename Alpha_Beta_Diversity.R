library(phyloseq)
library (ape)
library (ggplot2)

#----------------------------------- Alpha diversity analysis on complete data (Total-Data)
#-- Make a phyloseq object after filtering of OTUs
# Read in OTU table
otu_table_in <- read.csv("selected_feature-table_copy.txt", sep = "\t", row.names = 1)
total_asv <- colSums(otu_table_in)
range (total_asv)
sd (total_asv)
mean (total_asv)

# Read in taxonomy
# Separated by kingdom, phylum, class, order, family, genus, species
taxonomy <- read.csv("taxonomy.tsv", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)
# Read in metadata
metadata <- read.table("selected_Metadata.tsv", sep="\t", row.names = 1, header=T)
# Read in tree
phy_tree <- read_tree("tree.nwk")
# Import all as phyloseq objects
OTU <- otu_table(otu_table_in, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(metadata)
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

#Colors <- c("darkolivegreen4", "red", "darkgreen", "salmon4", "wheat", "goldenrod1", "mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine")
#Colors <- c("darkolivegreen4", "red", "darkgreen", "salmon4", "wheat", "goldenrod1", "mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine", "gray10", "gray15", "gray20", "gray25", "gray30", "gray35", "gray40", "gray45", "gray50", "gray55", "gray60")
#Colors <- c("wheat", "goldenrod1", "mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine", "gray10", "gray15", "gray20", "gray25", "gray30", "gray35", "gray40", "gray45", "gray50", "gray55", "gray60", "darkolivegreen4", "red", "darkgreen", "salmon4")
#Colors <- c("mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine")
#Colors <- c("wheat", "goldenrod1", "mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine")
Colors <- c("darkolivegreen4", "red", "salmon4", "darkgreen", "wheat", "goldenrod1")
#Colors1 <- c("wheat4", "goldenrod4", "mediumpurple4", "grey27", "blue", "paleturquoise4", "aquamarine4",  "gray10", "gray15", "gray20", "gray25", "gray30", "gray35", "gray40", "gray45", "gray50", "gray55", "gray60", "darkolivegreen","red4", "green", "salmon")
#Colors <- c("lightblue", "gray")

#Colors1 <- c("mediumpurple4", "grey27", "blue", "paleturquoise4", "aquamarine4")
Colors1 <- c("darkolivegreen", "red4", "salmon", "green", "wheat4", "goldenrod4")
#Colors1 <- c("wheat4", "goldenrod4", "mediumpurple4", "grey27", "blue", "paleturquoise4", "aquamarine4")
#Colors <- c("darkred", "darkgray")
#Colors1 <- c("midnightblue", "darkgray")
p <- plot_richness(ps, "diet", measures = c("Observed", "Chao1", "Shannon"), color = NULL, shape = NULL)
p <- p + geom_boxplot(data = p$data, aes(x= diet, y = value, fill = diet)) + 
  labs(x="",y="Alpha Diversity Measure") + 
  theme_bw() +
  scale_color_manual(values = Colors1)+
  scale_fill_manual(values = Colors)+
  geom_point(aes(colour = factor(diet)), position=position_jitterdodge(jitter.width = 1), size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 24)) 
ggsave("phyloseq_analysis-richness_estimates.pdf", p, width = 10, height = 10, limitsize = FALSE)
jpeg("Phyloseq_Richness.jpg", height = 10, width = 15, units = 'in', res = 600)
p
dev.off ()
erich <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
erich_diet <- cbind(erich, metadata)
write.table (erich_diet, file = "Richness_total.txt", sep = "\t")

kruskalmc(erich_diet$Observed ~ erich_diet$diet, probs=0.05)
wilcox.test(erich_diet$Shannon ~ erich_diet$diet, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)


#----------------------------To rarefy the OTU table
library("GUniFrac")
OTU_table <- read.table(file = "selected_features.txt", sep = "\t", header = T, row.names = 1)
otu.tab.rff <- Rarefy(OTU_table, depth = min(rowSums(OTU_table)))$otu.tab.rff # Samples in Row and OTUs in column
OTU_table_rar <- data.frame(otu.tab.rff)

write.table(OTU_table_rar, file = "Selected_FeatureTable_rarefied.txt", quote = FALSE, sep = '\t')


#---------------------------PCoA Analysis

Bushman2  = transform_sample_counts(ps, function(x) x / sum(x) )

#-------- weighted UniFrac PCoA
UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
UniFrac_distances <- as.matrix(UniFrac_distances)
UniFrac_dist_column <- melt(UniFrac_distances)
write.table (UniFrac_dist_column, file = "UniFrac_distances.txt", sep = "\t")

Uni_pcoa <- pcoa(UniFrac_distances)
Uni_pcoa$values[1:2,]
mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
pc <- c(1,2)

jpeg("UniFrac_PCoA.jpg", height = 10, width = 10, units = 'in', res = 600)
plot(Uni_pcoa$vectors[,1:2], bg= c("darkolivegreen4", "red", "salmon4")[as.factor(Bushman2@sam_data$diet)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$diet, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("darkolivegreen4", "red", "salmon4"))
ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$diet, lty=3, spider ="centroid", lwd=1, col="black")
#legend("bottomleft", legend = c("A_BL_Chow", "B_BL_WSD", "C_HF_Chow", "D_HF_WSD", "Gavage_V1", "Gavage_V3", "Selected_V1", "Selected_V3", "V1", "V3"), col = c("darkolivegreen4", "red", "darkgreen", "salmon4", "wheat", "goldenrod1", "mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
legend("bottomleft", legend = c("Chow", "HFD", "WSD"), col = c("darkolivegreen4", "red", "salmon4"), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
adonis2(UniFrac_distances ~ Bushman2@sam_data$diet)

#------------ unweighted UniFrac PCoA

#UniFrac_distances <- UniFrac(Bushman2, weighted=FALSE)
UniFrac_distances <- distance(Bushman2, method="bray")
UniFrac_distances <- as.matrix(UniFrac_distances)
UniFrac_dist_column <- melt(UniFrac_distances)
#write.table (UniFrac_dist_column, file = "unweighted_UniFrac_distances.txt", sep = "\t")
write.table (UniFrac_dist_column, file = "Bray_Curtis_distances.txt", sep = "\t")
Uni_pcoa <- pcoa(UniFrac_distances)
Uni_pcoa$values[1:2,]
mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
pc <- c(1,2)

#jpeg("unweighted_UniFrac_PCoA.jpg", height = 10, width = 10, units = 'in', res = 600)
jpeg("Bray_ Curtis_PCoA.jpg", height = 10, width = 10, units = 'in', res = 600)
plot(Uni_pcoa$vectors[,1:2], bg= c( "mediumpurple3", "snow4", "deepskyblue", "aquamarine")[as.factor(Bushman2@sam_data$diet)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$diet, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("mediumpurple3", "snow4", "deepskyblue", "aquamarine"))
#ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$diet, lty=3, spider ="centroid", lwd=1, col="black")
legend("topleft", legend = c("Selected_V1", "Selected_V3", "V1", "V3"), col = c("mediumpurple3", "snow4", "deepskyblue", "aquamarine"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
#legend("bottomleft", legend = c("A20", "control", "HHSA", "HSA", "LSA", "SA", "Thr"), pch=c(21,22),cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
adonis2(UniFrac_distances ~ Bushman2@sam_data$diet)
