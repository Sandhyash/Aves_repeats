a=read.table("simp_d_orderwise.txt",header=F)
library(stringr)
a$V1=str_to_title(a$V1)
jpeg("simpson_diversity_aves_repeats.jpeg",width=23.2,height=13,units="in",res=300)
par(mai=c(2.6,2.6,0.5,0.2),mgp=c(11,0.2,0))
barplot(a$V2,names=a$V1,ylab="",ylim=c(0,1),xlab="",col="#F7B37E",main="Simpson diversity of repeats across different orders",las=3,cex.axis=1,cex.lab=1.5)
title(xlab="Orders",line=8,cex.lab=2)
title(ylab="Simpson diversity",cex.lab=2,line=2)
dev.off()
q()
n

