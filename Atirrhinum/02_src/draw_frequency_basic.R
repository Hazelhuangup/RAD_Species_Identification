args<-commandArgs(TRUE)

pdf(args[1],width=8,height=6)
par(mar=c(4,2,4,2))
file <- read.table(args[2],sep="\t")
plot(file$V2,file$V1,type='h',pch=20,xlab="No_individuals_not_missing",ylab="No_loci",lwd=2)
dev.off()
q()
