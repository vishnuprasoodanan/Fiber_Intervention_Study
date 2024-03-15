library(phyloseq)
library(ape)
library(vegan)
library(reshape2)
library(ggplot2)

#----------------------------------- Alpha diversity analysis on complete data (Total-Data)
#-- Make a phyloseq object after filtering of OTUs
# Read in OTU table
otu_table_in <- read.csv("selected_feature-table_copy.txt", sep = "\t", row.names = 1)

# Read in taxonomy
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

# Same sample names
sample_names(OTU) <- gsub('[.]', '-', sample_names(OTU))
sample_names(OTU) <- gsub('^[X]', '', sample_names(OTU)) 
head(sample_names(OTU))
head(sample_names(META))

# create Phyloseq object
ps <- phyloseq(OTU, TAX, META, phy_tree)

Colors <- c("darkolivegreen4", "red", "darkgreen", "salmon4")
Colors1 <- c("darkolivegreen", "red4", "green", "salmon")

p <- plot_richness(ps, "diet", measures = c("Observed", "Chao1", "Shannon"), color = NULL, shape = NULL)
p <- p + geom_boxplot(data = p$data, aes(x= diet, y = value, fill = diet)) + 
  labs(x="",y="Alpha Diversity Measure") + 
  theme_bw() +
  scale_color_manual(values = Colors1)+
  scale_fill_manual(values = Colors)+
  geom_point(aes(colour = factor(diet)), position=position_jitterdodge(jitter.width = 1), size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 24)) 
ggsave("phyloseq_analysis-richness_estimates.pdf", p, width = 10, height = 10, limitsize = FALSE)
dev.off ()
erich <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
erich_status <- cbind(erich, metadata)
write.table (erich_status, file = "Richness_total.txt", sep = "\t")

kruskalmc(erich_staus$Observed ~ erich_status$status, probs=0.05)
wilcox.test(erich_status$Shannon ~ erich_status$status, paired = FALSE, p.adjust.method="fdr", alternative = c("two.sided"), conf.level = 0.95)


#---------------------------PCoA Analysis

otu_rel_ab  = transform_sample_counts(ps, function(x) x / sum(x) )

#-------- weighted UniFrac PCoA
UniFrac_distances <- UniFrac(otu_rel_ab, weighted=TRUE)
UniFrac_distances <- as.matrix(UniFrac_distances)
UniFrac_dist_column <- melt(UniFrac_distances)
write.table (UniFrac_dist_column, file = "UniFrac_distances.txt", sep = "\t")

Uni_pcoa <- pcoa(UniFrac_distances)
Uni_pcoa$values[1:2,]
mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
pc <- c(1,2)

jpeg("UniFrac_PCoA.jpg", height = 10, width = 10, units = 'in', res = 600)
plot(Uni_pcoa$vectors[,1:2], bg= c("darkolivegreen4", "red", "darkgreen", "salmon4")[as.factor(otu_rel_ab@sam_data$status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
text(Uni_pcoa$vectors[,1:2], labels=rownames(otu_rel_ab@sam_data), cex=0.3, font=1, pos=1)
ordiellipse(Uni_pcoa$vectors[,1:2], otu_rel_ab@sam_data$status, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("darkolivegreen4", "red", "darkgreen", "salmon4"))
ordispider(Uni_pcoa$vectors[,1:2], otu_rel_ab@sam_data$status, lty=3, spider ="centroid", lwd=1, col="black")
legend("bottomleft", legend = c("A_BL_Chow",  "B_BL_WSD", "C_HF_Chow", "D_HF_WSD"), col = c("darkolivegreen4", "red", "darkgreen", "salmon4"), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
adonis2(UniFrac_distances ~ otu_rel_ab@sam_data$status)

#------------ unweighted UniFrac PCoA and Bray-Curtis distance calculations
#calculate unweighted unifrac distance
#UniFrac_distances <- UniFrac(otu_rel_ab, weighted=FALSE)

#calculate Bray-Curtis distance
UniFrac_distances <- distance(otu_rel_ab, method="bray")
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
plot(Uni_pcoa$vectors[,1:2], bg= c( "mediumpurple3", "snow4", "deepskyblue", "aquamarine")[as.factor(otu_rel_ab@sam_data$status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#text(Uni_pcoa$vectors[,1:2], labels=rownames(otu_rel_ab@sam_data), cex=0.3, font=1, pos=1)
ordiellipse(Uni_pcoa$vectors[,1:2], otu_rel_ab@sam_data$status, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("mediumpurple3", "snow4", "deepskyblue", "aquamarine"))
#ordispider(Uni_pcoa$vectors[,1:2], otu_rel_ab@sam_data$status, lty=3, spider ="centroid", lwd=1, col="black")
legend("topleft", legend = c("Selected_V1", "Selected_V3", "V1", "V3"), col = c("mediumpurple3", "snow4", "deepskyblue", "aquamarine"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
#legend("bottomleft", legend = c("A20", "control", "HHSA", "HSA", "LSA", "SA", "Thr"), pch=c(21,22),cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
abline(h=0, v=0, col = "gray60")
dev.off ()
adonis2(UniFrac_distances ~ otu_rel_ab@sam_data$status)
