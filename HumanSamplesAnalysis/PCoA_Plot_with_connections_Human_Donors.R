#

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
plot(Uni_pcoa$vectors[,1:2], bg= c("white", "white", "pink1", "purple1", "goldenrod1", "grey75", "steelblue")[as.factor(Bushman2@sam_data$Status1)], pch=21, cex=2, xlab=paste0("PCoA", pc[1], " (", mds.var.per[1], "%)"), ylab=paste0("PCoA", pc[2], " (", mds.var.per[2], "%)"))
ordiellipse(Uni_pcoa$vectors[,1:2], Bushman2@sam_data$Status, kind="sd", lwd=1, lty=3, draw = "polygon", alpha = 60, col = c("skyblue","palegreen", "orangered"))

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
