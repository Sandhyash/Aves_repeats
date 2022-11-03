### Frequency barplots for repeat frequency in aves clade
b<-read.table("aves_morethanseventy_combined")
b1<-data.frame(table(b$V11))
b2<-b1[order(b1$Freq,decreasing=T),]
write.table(b2,file="repeat_frequency_aves.txt",quote=F,row.names=F,col.names=F)
b3<-b2[1:30,]
pdf("Figure_S4.pdf")
barplot(b3$Freq,names=b3$Var1,las=2,col="plum",xlab="Amino acid repeats",ylab="Frequency",ylim=c(0,10000),main="Frequency of amino acid repeats in Aves")
dev.off()
jpeg("Figure_S4.jpeg",width=28,height=18,res=300,units="in")
par(mai=c(2,2,1.5,1.5))
barplot(b3$Freq,names=b3$Var1,las=2,col="plum",xlab="",ylab="",ylim=c(0,10000),main="",cex.axis=2,cex.names=2)
title(main="Frequency of amino acid repeats in Aves" ,line=1,cex.main=2.3)
title(ylab="Frequency", line=4, cex.lab=2.5)        
title(xlab="Amino acid repeats", line=5, cex.lab=2.5)
dev.off()

