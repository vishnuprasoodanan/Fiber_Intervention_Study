#Plotting the boxplot of relative abundance/clr-transformd count data of each group (HD-WSD, HD-Chow, HF-WAS, and HF-Chow)
#Input file contain the abundance of each sample. Taxa (Genus or Phylum) in the rows and samples in the colums. The first coulumn (named 'Status') contain group information of each sample

data <- read.table(file = "boxplot_input.txt", sep = "\t", header = T)
pdf(file = "Abundance_boxplots.pdf")
colnames(data) -> Species_name


A_BL_Chow <- data[data$Status == "AK_BL_Chow_termination",]
B_BL_WSD <- data[data$Status == "AL_BL_WSD_termination",]
C_HF_Chow <- data[data$Status == "AM_HF_Chow_termination",]
C_HF_WSD <- data[data$Status == "AN_HF_WSD_termination",]


Colors <- c("darkolivegreen4", "red", "darkgreen", "salmon4")
Colors1 <- c("darkolivegreen","red4", "green", "salmon")

plot_lst <- vector("list", length = 23)
for (i in 3:ncol(data)) {    
  
  species = Species_name[i]
  data1 = data[,c(1:2,i)]
  colnames(data1) <- c("ID", "Status", "name")
  
  A_BL_Chow_all <- A_BL_Chow[,i]
  B_BL_WSD_all <- B_BL_WSD[,i]
  C_HF_Chow_all <- C_HF_Chow[,i]
  C_HF_WSD_all <- C_HF_WSD[,i]
  
  
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

#Kruskal-Wallis test for the genera
sel_genera <- c("Blautia", "Bacteroides", "Akkermansia", "Clostridium_innocuum_group", "Lachnoclostridium", "Bilophila", "Coprobacillus", "Enterococcus")
for (i in sel_genera) {
  # Create the formula for the test
  formula <- paste(i, "~ Status")
  
  # Run the kruskalmc test
  result <- kruskalmc(as.formula(formula), data = data, probs = 0.01, alpha = 0.01)
  
  # Print the result
  print(paste("test for",i))
  print(result)
}
