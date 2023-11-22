# Fiber_Intervention_Study

### ASV data used for analysis
**SampleWise_ASV_Summary.xlsx:** Samplewise summary of ASV table contaning total number of ASVs detected per sample and the number of distint ASVs detected per sample

**ASV_ID_Summary.xlsx:** Number of sample each ASVs present and total number of ASVs detected in all samples



**Pairwise_Wilcox_Test.R:** This script is used to carry out pairwise Wilcoxon test in genus abundance data from visit-1 and visit-3.

**PCoA_Plot_with_connections_Human_Donors.R:** This script is used to plot PCoA using any intersample distance matrix (weighted/unweighted unifrac, Bray-Curtis). This script also connect the paired samples (samples from visit-1 and visit-3) using solid lines. This paired information should be mentioned in metadata file (Status2 column).

**Box_Plot_Genera_Abundance.R:** This script is used to plot box-plots based on the relative abundance information of microbial genera identified in a dataset. Sample-names will be in rows and genera will be columns of the input file. The last column named "Status" contain grouping information

**Stacked_BarPlot_Genus_Abundance.R:**  First part of this script is used to plot stacked bar-plots based on the relative abundance information of microbial genera identified in a dataset. Second part of this script will plot the stacked bar-plot by calculating the average relative abundance of genera in each group of sample.

**Relative_Abundance.R**



