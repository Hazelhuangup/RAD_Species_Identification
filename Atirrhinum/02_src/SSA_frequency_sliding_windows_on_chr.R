#!/usr/bin/Rscript
args<-commandArgs(TRUE)
dir <- args[1]
setwd(dir)
freq_file_list <- as.vector(Sys.glob("*freq"))

# Libraries
library(ggplot2)

pdf('Atirrhinum.SSA.sliding_win.freq_all.pdf',width=40,height=6, onefile = TRUE)
for(freq_file in freq_file_list){
	data <- read.table(freq_file, NULL, header = FALSE)
	mydf <- data.frame(
		Chr_Name= data$V1,
		Position = as.numeric(data$V2),
		Frequency = as.numeric(data$V4)
)
	plot <-
	ggplot(mydf, aes(fill = Chr_Name, col = Chr_Name, x = Position, y = Frequency), group = Chr_Name) + 
		geom_point() +
		ggtitle(as.character(freq_file)) +
		theme(plot.title = element_text(hjust = 0.5)) +
		facet_grid( ~Chr_Name) +
		theme_light()
	print(plot)
#	ggsave(paste(freq_file,"pic.pdf",sep=''), width = 40, height = 6)
}
dev.off()
