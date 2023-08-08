library(phyloseq)
library (ape)
library (ggplot2)
#----------------------------------- Alpha diversity analysis on complete data (Total-Data)
#-- Make a phyloseq object after filtering of OTUs
# Read in OTU table
otu_table_in <- read.csv("selected-feature-table_EYFD.tsv", sep = "\t", row.names = 1)
total_asv <- colSums(otu_table_in)
range (total_asv)
sd (total_asv)
mean (total_asv)
# Read in taxonomy
# Separated by kingdom, phylum, class, order, family, genus, species
taxonomy <- read.csv("taxonomy.tsv", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)
# Read in metadata
metadata <- read.table("Metadata_V2.tsv", sep="\t", row.names = 1, header=T)
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
Bushman2  = transform_sample_counts(ps, function(x) x / sum(x) )
#-------- weighted UniFrac PCoA
UniFrac_distances <- UniFrac(Bushman2, weighted=TRUE)
#UniFrac_distances <- distance(Bushman2, method="bray", type="samples")
UniFrac_distances <- as.matrix(UniFrac_distances)
UniFrac_dist_column <- melt(UniFrac_distances)
write.table (UniFrac_dist_column, file = "Unweighted_UniFrac_distances.txt", sep = "\t")
Uni_pcoa <- pcoa(UniFrac_distances)
Uni_pcoa$values[1:2,]
mds.var.per = round(Uni_pcoa$values$Eigenvalues/sum(Uni_pcoa$values$Eigenvalues)*100, 1)
pc <- c(1,2)
jpeg("BrayCurtis_check_v2.jpg", height = 10, width = 10, units = 'in', res = 600)
pcoa_df <- Uni_pcoa$vectors[,1:2]
#plot(Uni_pcoa$vectors[,1:2], bg= c("darkolivegreen4", "red", "darkgreen", "salmon4", "darkolivegreen", "red3","forestgreen" , "salmon")[as.factor(Bushman2@sam_data$Status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("white", "white", "pink1", "purple1", "goldenrod1", "grey75", "orangered", "steelblue")[as.factor(Bushman2@sam_data$Status1)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("white", "white", "pink1", "goldenrod1", "grey75", "orangered", "steelblue")[as.factor(Bushman2@sam_data$Status1)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("pink1", "goldenrod1", "grey75", "orangered", "steelblue")[as.factor(Bushman2@sam_data$Status1)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("pink1", "purple1", "goldenrod1", "grey75", "steelblue")[as.factor(Bushman2@sam_data$Status1)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
plot(Uni_pcoa$vectors[,1:2], bg= c("white", "white", "pink1", "purple1", "goldenrod1", "grey75", "steelblue")[as.factor(Bushman2@sam_data$Status1)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("pink1", "purple1", "goldenrod1", "grey75", "orangered", "steelblue")[as.factor(Bushman2@sam_data$Status1)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("steelblue", "sienna1", "deepskyblue")[as.factor(Bushman2@sam_data$Status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("red","darkolivegreen4", "orange", "darkgreen", "salmon4", "limegreen")[as.factor(Bushman2@sam_data$Status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("orange", "darkgreen", "salmon4", "limegreen")[as.factor(Bushman2@sam_data$Status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("Blue", "PowderBlue", "orangered", "red4")[as.factor(Bushman2@sam_data$Status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("cyan", "midnightblue", "chocolate", "salmon4")[as.factor(Bushman2@sam_data$Status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#plot(Uni_pcoa$vectors[,1:2], bg= c("darkolivegreen4", "red", "darkgreen", "salmon4")[as.factor(Bushman2@sam_data$Status)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
#text(Uni_pcoa$vectors[,1:2], labels=rownames(Bushman2@sam_data), cex=0.3, font=1, pos=1)
#ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$Status2, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 50, col = c("skyblue","palegreen"))
ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$Status, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 60, col = c("skyblue","palegreen", "orangered"))
#ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$Status, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 60, col = c("steelblue", "sienna1", "deepskyblue"))
#ordispider(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$Status, lty=3, spider ="centroid", lwd=1, col="black")
#legend("topright", legend = c("AK_BL_Chow_termination", "AL_BL_WSD_termination", "AM_HF_Chow_termination", "AN_HF_WSD_termination"), col = c("darkolivegreen4", "red", "darkgreen", "salmon4", "wheat", "goldenrod1", "mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
#legend("bottomleft", legend = c("AQ_BL_Habitual_diet_Citro", "AR_BL_High_fiber_Citro"), col = c("darkolivegreen4", "red", "darkgreen", "salmon4", "wheat", "goldenrod1", "mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
#legend("bottomleft", legend = c("BL_Bcoccoides", "BL_Bcoccoides_control", "T_Bcoccoides", "T_Bcoccoides_control"), col = c("cyan", "midnightblue", "chocolate", "salmon4"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")
#legend("bottomright", legend = c("Baseline", "FMT-Termination", "Human"), col = c("steelblue", "sienna1", "deepskyblue"),lty = c(1,1,1,1), cex=0.7, title = "", border = "white", fill = NULL, bg="white", bty = "n")

for (i in 1:(length(Bushman2@sam_data$Status2)-1)) {
if (!is.na(Bushman2@sam_data$Status2[i])) {
   for (j in (i+1):length(Bushman2@sam_data$Status2)) {
     if (!is.na(Bushman2@sam_data$Status2[j]) && Bushman2@sam_data$Status2[i] == Bushman2@sam_data$Status2[j]) {
       x <- c(pcoa_df[i, "Axis.1"], pcoa_df[j, "Axis.1"])
       y <- c(pcoa_df[i, "Axis.2"], pcoa_df[j, "Axis.2"])
       lines(x, y, col = "black", lty = "solid")
      }
    }
  }
}

abline(h=0, v=0, col = "gray60")
dev.off ()
