##-------------------------------------------------------------------------------
#boxplot for core microbiome genera
#location: ~/WORK/FIBER_INTERVENTION_16S/CORE_MICROBIOME_STACKED_BARPLOTS/BOXPLOTS
data <- read.table(file = "Group_data.txt", sep = "\t", header = T)
pdf(file = "25core_microbiome_Rel_Abundance_boxplots.pdf")
colnames(data) -> Species_name

A_BL_Chow <- data[data$Status == "A_BL_Chow",]
B_BL_WSD <- data[data$Status == "B_BL_WSD",]
C_HF_Chow <- data[data$Status == "C_HF_Chow",]
C_HF_WSD <- data[data$Status == "C_HF_WSD",]
Gavage_V1 <- data[data$Status == "Gavage_V1",]
Gavage_V3 <- data[data$Status == "Gavage_V3",]
Selected_V1 <- data[data$Status == "Selected_V1",]
Selected_V3 <- data[data$Status == "Selected_V3",]
V1 <- data[data$Status == "V1",]
V3 <- data[data$Status == "V3",]

Colors <- c("darkolivegreen4", "red", "darkgreen", "salmon4", "wheat", "goldenrod1", "mediumpurple3", "snow4", "deepskyblue", "paleturquoise" ,"aquamarine")
Colors1 <- c("darkolivegreen","red4", "green", "salmon", "wheat4", "goldenrod4", "mediumpurple4", "grey27", "blue", "paleturquoise4", "aquamarine4")

plot_lst <- vector("list", length = 23)
for (i in 3:ncol(data)) {    
  
  species = Species_name[i]
  data1 = data[,c(1:2,i)]
  colnames(data1) <- c("ID", "Status", "name")
  
  A_BL_Chow_all <- A_BL_Chow[,i]
  B_BL_WSD_all <- B_BL_WSD[,i]
  C_HF_Chow_all <- C_HF_Chow[,i]
  C_HF_WSD_all <- C_HF_WSD[,i]
  Gavage_V1_all <- Gavage_V1[,i]
  Gavage_V3_all <- Gavage_V3[,i]
  Selected_V1_all <- Selected_V1[,i]
  Selected_V3_all <- Selected_V3[,i]
  V1_all <- V1[,i]
  V3_all <- V3[,i]
  
  print(species)
  species <- gsub("_", " ", species)
  P<- ggplot(data1, aes(Status, name, fill=Status))+
    ggtitle(species)+
    labs(y = "Relative-Abundance")+
    geom_boxplot(outlier.shape=NA)+ 
    geom_point(aes(colour = factor(data$Status)), position=position_jitterdodge(jitter.width = 0.5))+
    scale_color_manual(values = Colors1)+
    scale_fill_manual(values = Colors)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot_lst[[i]] <- P
}
ml <- marrangeGrob(plot_lst, nrow = 1, ncol = 1)
print(ml)
dev.off()
