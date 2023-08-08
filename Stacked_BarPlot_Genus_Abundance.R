library(ggplot2)
df <- read.table(file="Mice.txt", sep = "\t", header = T)
data_order <- df[order(apply(df, 1, min)),]
# convert data from wide to long format using tidyr
library(tidyr)
df_long <- pivot_longer(df, cols = -Name, names_to = "Sample", values_to = "Abundance")
#colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "deepskyblue", "goldenrod1", "darkolivegreen", "salmon", "mediumpurple3", "paleturquoise4", "goldenrod4", "#FF00FF", "#F433FF", "#E600FF", "#CC00FF", "#B266FF","#9B30FF", "#8800FF", "#660066", "#993399", "#CC66CC", "#FF33CC", "#FF0099", "gray")
#colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "deepskyblue", "goldenrod1", "darkolivegreen", "salmon", "mediumpurple3", "paleturquoise4", "goldenrod4", "gray")
#colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "deepskyblue", "gray")
colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "gray")
# create stacked barplot using ggplot2
ggplot(df_long, aes(x = Sample, y = Abundance, fill = Name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "Abundance of Features in Each Sample", x = "Sample", y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

###--------------------------------------------------------------------------
### plot stacked bar-plots using the average relative abundance
library(ggplot2)
library(dplyr)
set.seed(123)
data <- read.table(file="group_data.txt", sep = "\t", header = T)
# Compute average abundance by group
# group_averages <- data %>% group_by(Group) %>% 
#   summarise(across(c("B_Bacteroides", "E_Lachnospiraceae_unknowngenus", "A_Alistipes", "C_Blautia", "D_Lachnoclostridium", "F_Butyricicoccus", "G_Christensenellaceae_R.7_group", "H_Incertae_Sedis", "I_Bilophila", "J_Oscillibacter", "K_Phascolarctobacterium", "L_Parabacteroides", "M_Coprococcus", "N_Subdoligranulum", "O_Anaerostipes", "P_Bifidobacterium", "Q_Colidextribacter", "R_Barnesiella", "S_Parasutterella", "T_Akkermansia", "U_Gastranaerophilales", "V_Ruminococcus_torques_group", "W_Sellimonas", "X_Other"), ~ mean(.x, na.rm = TRUE)))
group_averages <- data %>% group_by(Group) %>% 
  summarise(across(c("B_Bacteroides", "E_Lachnospiraceae_unknowngenus", "A_Alistipes", "C_Blautia", "D_Lachnoclostridium", "F_Butyricicoccus", "G_Christensenellaceae_R.7_group", "H_Incertae_Sedis", "I_Bilophila", "J_Oscillibacter", "K_Phascolarctobacterium", "L_Other"), ~ mean(.x, na.rm = TRUE)))

write.table(group_averages, file = "group_average_rel_abundance.txt", sep = '\t', quote = FALSE, row.names = FALSE)

df <- read.table(file="group_average_rel_abundance.txt", sep = "\t", header = T)
df_long <- pivot_longer(df, cols = -Group, names_to = "Sample", values_to = "Abundance")
#colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#C2C2C2",
           # "#17becf", "#1f77b4", "#9467bd", "#d62728", "#2ca02c", "#ff7f0e", "#8c564b", "#bcbd22", "#7f7f7f",
           # "#dbdb8d", "#ffbb78", "#98df8a", "#ff9896", "#f7b6d2", "#c5b0d5", "#8c8c8c")
#colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "deepskyblue","gray")
colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "deepskyblue", "goldenrod1", "darkolivegreen", "salmon", "mediumpurple3", "paleturquoise4", "goldenrod4", "gray")
#colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "deepskyblue", "goldenrod1", "darkolivegreen", "salmon", "mediumpurple3", "paleturquoise4", "goldenrod4", "#FF00FF", "#F433FF", "#E600FF", "#CC00FF", "#B266FF","#9B30FF", "#8800FF", "#660066", "#993399", "#CC66CC", "#FF33CC", "#FF0099", "gray")
# create stacked barplot using ggplot2
ggplot(df_long, aes(x = Group, y = Abundance, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors) +
  labs(title = "Abundance of Features in Each Sample", x = "Group", y = "Average Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
