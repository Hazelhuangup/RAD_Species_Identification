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
ggplot(mydf, group = Seq_Method, aes(fill = Seq_method, col = Seq_method, xmin = x, xmax = x+0.6, ymin = 0, ymax = Monophyly_Ratio)) + 
  geom_rect() +
  theme_light()

