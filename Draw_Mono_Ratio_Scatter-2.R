# Working directory
setwd('/Users/Mac/Documents/PhD Research/Drawing')

# Libraries
library(ggplot2)
data <- read.csv('RAD_published_species.txt', NULL, header = TRUE)
mydf <- data.frame(
  Specie_Name= data$Material, 
  Monophyly_Ratio = as.numeric(data$Monophyletic_Ratio),
  TaxonCategory = data$Kingdom,
  NoSpecies = as.numeric(data$No.N.1),
  Seq_method = data$Seq_Method,
  x = data$x
)
ggplot(mydf, aes(fill = TaxonCategory, col = TaxonCategory, x = NoSpecies, y = Monophyly_Ratio)) + 
  geom_jitter(size=3) +
  theme_light()
ggplot(mydf, aes(fill = Seq_method, col = Seq_method, x = x , y = Monophyly_Ratio), group = Seq_method) + 
  geom_point(size=3) +
  geom_line(size=1) +
  theme_light()
