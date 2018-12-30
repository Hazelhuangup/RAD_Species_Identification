# Working directory
setwd('~/Desktop/PhD Research/Drawing/')

# Libraries
library('ggplot2')
data <- read.csv('RAD_published_species.txt', NULL, header = TRUE)
mydf <- data.frame(
  Specie_Name= head(data$Material, n=39L), 
  Monophyly_Ratio = as.numeric(head(data$Monophyletic_Ratio,n=39L)),
  TaxonCategory = head(data$Kingdom,n=39L),
  NoSpecies = as.numeric(head(data$No.N.1, n=39L))
)
ggplot(mydf, aes(col = TaxonCategory, x = NoSpecies, y = Monophyly_Ratio)) + 
  geom_point(size=3) +
  theme_light()

