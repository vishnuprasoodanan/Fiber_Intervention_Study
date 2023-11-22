# Fiber_Intervention_Study
### Analysis Summary

The analysis involved sequencing data from 142 human participants, 207 mouse samples associated with different experiments and three control samples. A total of 19,565,524 (median = 91,452) paired-end reads from mouse samples and 4,764,753 (median = 32,206) paired-end reads from human samples were used for the analyses (Sheet-1,2, Supplementary Table-1). The reads were subjected to quality filtration using Fastx (http://hannonlab.cshl.edu/fastx_toolkit/index.html) by the criteria that sequences with at least 70% of bases with the quality score of ≥ 25. Assessed sequence quality and adapter presence using FastQC (Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
), and then proceeded to remove the adapters. Analysis of the human and mouse-derived 16S amplicon sequences was performed using the QIIME21 (version 2021.4) pipeline and R (version 4.1.3) in R Studio (RStudio Team, version 2022.07.2). Sequences were then clustered into ASVs using the SILVA2 classifier (version 138). DADA2 was employed to generate the feature table, which comprises ASV-IDs along with the corresponding count of detected ASVs in each sample. 9,974 ASVs were detected in this analysis. Low-abundance ASVs were excluded based on the criterion that an ASV should be present in a minimum of two samples with a total ASV count of less than 10. 3,115 ASVs remained after this criterion. On average 23,900 ASVs were detected per samples (median =13,183).

### ASV data generated from QIIME2 analysis
The ASV table generated from this analysis has been stored in the directory named "ASV_Table_Summary". Within this folder, you'll also find a README file that provides details about the contents of each file included.

### Scripts used for Analysis
**Pairwise_Wilcox_Test.R:** This script is used to carry out pairwise Wilcoxon test in genus abundance data from visit-1 and visit-3.

**PCoA_Plot_with_connections_Human_Donors.R:** This script is used to plot PCoA using any intersample distance matrix (weighted/unweighted unifrac, Bray-Curtis). This script also connect the paired samples (samples from visit-1 and visit-3) using solid lines. This paired information should be mentioned in metadata file (Status2 column).

**Box_Plot_Genera_Abundance.R:** This script is used to plot box-plots based on the relative abundance information of microbial genera identified in a dataset. Sample-names will be in rows and genera will be columns of the input file. The last column named "Status" contain grouping information

**Stacked_BarPlot_Genus_Abundance.R:**  First part of this script is used to plot stacked bar-plots based on the relative abundance information of microbial genera identified in a dataset. Second part of this script will plot the stacked bar-plot by calculating the average relative abundance of genera in each group of sample.

**Relative_Abundance.R**



