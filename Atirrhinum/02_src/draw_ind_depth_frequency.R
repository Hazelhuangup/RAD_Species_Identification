#!/usr/bin/Rscript
args<-commandArgs(TRUE)
dir <- args[1]
setwd(dir)
freq_file_list <- as.vector(Sys.glob("*freq"))

# Libraries
library(ggplot2)

pdf('Atirrhinum.read_depth.freq_all.pdf',width=10,height=6, onefile = TRUE)
for(freq_file in freq_file_list){
            data <- read.table(freq_file, NULL, header = FALSE)
            mydf <- data.frame(
                        depth = data$V2,
                        Frequency = as.numeric(data$V1)
)
            plot <-
            ggplot(mydf, aes(col = depth, x = depth, y = Frequency)) +
						scale_y_continuous(trans='log10') +
                        geom_point() +
                        ggtitle(as.character(freq_file)) +
                        theme(plot.title = element_text(hjust = 1)) +
                        theme_light()
            print(plot)
}
dev.off()
