#Script for generating CLR-transformed abundance values from Genus (or Phylum) abundance/count table

#"FMT_Genus_count_table.txt" contain the number of ASVs assigned to each genera/phylum
otu_table_in <- read.csv("FMT_Genus_count_table.txt", sep = "\t", row.names = 1)
# Read in taxonomy file containing Genus ID and complete genus annotations
taxonomy <- read.csv("FMT_Taxonomy.txt", sep = "\t", row.names = 1)
taxonomy <- as.matrix(taxonomy)
# Read in metadata
metadata <- read.table("FMT_metadata.txt", sep="\t", row.names = 1, header=T)

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
ps <- phyloseq(OTU, TAX, META)

genus_clr <- tax_transform(ps, trans = "clr", rank = "Taxon")

max_taxa <- 100
best_taxa <- names(sort(taxa_sums(genus_clr@otu_table), decreasing=TRUE)[1:max_taxa])

# get best max_taxa from the centered log ratio object
phylo_thresholded <- prune_taxa(best_taxa, genus_clr@otu_table)
write.table(phylo_thresholded, file = "clr-transformed_genus_abundance_100.txt", sep = '\t', quote = FALSE)
