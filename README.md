# Fiber_Intervention_Study
### Analysis Summary
The analysis involved sequencing data from 142 human participants, 207 mouse samples associated with different experiments and three control samples. A total of 19,536,751(median = 91,440) paired-end reads from mouse samples and 4,764,753 (median = 32,206) paired-end reads from human samples were used for the analyses (Sheet-1,2, Supplementary Table-1). The reads were subjected to quality filtration using Fastx (http://hannonlab.cshl.edu/fastx_toolkit/index.html) by the criteria that sequences with at least 70% of bases with the quality score of ≥ 25. Assessed sequence quality and adapter presence using FastQC (Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
), and then proceeded to remove the adapters. Analysis of the human and mouse-derived 16S amplicon sequences was performed using the QIIME21 (version 2021.4) pipeline and R (version 4.1.3) in R Studio (RStudio Team, version 2022.07.2). Sequences were then clustered into ASVs using the SILVA2 classifier (version 138). DADA23 was employed to generate the feature table, which comprises ASV-IDs along with the corresponding count of detected ASVs in each sample. 9,975 ASVs were detected in this analysis. Low-abundance ASVs were excluded based on the criterion that an ASV should be present in a minimum of two samples with a total ASV count of less than 10. 3,116 ASVs remained after this criterion. On average ~23,887 ASVs were detected per samples (median =13,183). Alpha and beta diversity analysis were carried out using the abundance of 3205 ASVs. The rarefaction analysis was performed using GUniFrac (https://cran.r-project.org/web/packages/GUniFrac/GUniFrac.pdf) package in R. The α-diversity metrics (observed species, Shannon, Simpson, and Chao1) and β-diversity (unweighted UniFrac distance, weighted UniFrac distance, and Bray–Curtis distance) were calculated using phyloseq4, vegan (https://cran.r-project.org/web/packages/vegan/vegan.pdf) and ape (https://cran.r-project.org/web/packages/ape/ape.pdf). The 'mia', 'miaViz', and 'vegan' packages were employed to compute the centered-log ratio transformation of ASV abundance data and assess intersample distances utilizing the Euclidean distance metric.
The number of reads assigned to different taxonomic classes (mainly phylum, genus, and species) was calculated, and the taxonomic composition was evaluated for each sample. The ASV count per taxonomic clade per sample was normalized by dividing it with the total number of reads in the corresponding sample. taxa_transform function in phyloseq was used to transform genus and phylum abundance data into clr-transformed values. Statistical differences between groups were calculated by using Kruskal–Wallis or PERMANOVA and 999 permutations (beta-diversity). The core microbiome (at the genus level)5 was determined based on the presence of genera in 25%, 50%, 75%, and 80% of samples within each group (human participants in visit-1 and visit-3, gavage, and mice received FMT). The plots were generated using ggplot26.


### ASV data generated from QIIME2 analysis
The ASV table generated from this analysis has been stored in the directory named "ASV_Table_Summary". Within this folder, you'll also find a README file that provides details about the contents of each file included.

### Scripts used for Analysis
**Pairwise_Wilcox_Test.R:** This script is used to carry out pairwise Wilcoxon test in genus abundance data from visit-1 and visit-3.

**PCoA_Plot_with_connections_Human_Donors.R:** This script is used to plot PCoA using any intersample distance matrix (weighted/unweighted unifrac, Bray-Curtis). This script also connect the paired samples (samples from visit-1 and visit-3) using solid lines. This paired information should be mentioned in metadata file (Status2 column).

**Box_Plot_Genera_Abundance.R:** This script is used to plot box-plots based on the relative abundance information of microbial genera identified in a dataset. Sample-names will be in rows and genera will be columns of the input file. The last column named "Status" contain grouping information

**Stacked_BarPlot_Genus_Abundance.R:**  First part of this script is used to plot stacked bar-plots based on the relative abundance information of microbial genera identified in a dataset. Second part of this script will plot the stacked bar-plot by calculating the average relative abundance of genera in each group of sample.

**Relative_Abundance.R**



