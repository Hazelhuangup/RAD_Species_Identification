# This script is for visualize the pattern of monophyletic ratio of species groups
# HazelHuang <W.Huang-24@sms.ed.ac.uk>

# Working directory
setwd('~/Desktop/PhD Research/Drawing/')

# 
data <- read.csv('Comprehensive information-Published groups with RAD-seq.csv')
TaxonCategory <- data$Kingdom
mydf <- data.frame(Specie_Name= data$Material, Monophyly_Ratio = data$Monophyletic_Ratio)
ggplot(mydf, aes(Specie_Name, Monophyly_Ratio)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
