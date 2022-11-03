a<-read.table("top10_rpt_prop_combined",header=T)
a$Repeats<-as.character(a$Repeats)
a[14,1]<-"Remaining" #add this row at the end
for(i in 2:length(a)){
a[14,i]<-1-sum(na.omit(a[,i]))
}
a<-data.frame(a,row.names=1)
cols=c("yellowgreen","violetred","violet","turquoise1","tomato","thistle","tan","steelblue1","springgreen3","slateblue2","salmon","purple","plum1","saddlebrown")
pdf("top10_repeat_proportion.pdf",width=30,height=20)
par(mai=c(3,3,1,3),mgp=c(13,0.9,0))
barplot(as.matrix(a),col=cols,main="Proportion of top 10 repeats",ylab="Order",xlab="Proportion of repeats",border=F,cex.axis=2,cex.names=1.5,cex.lab=1.5,horiz=T,las=1)
legend("topright",legend=rownames(a),fill=cols,bty="n",xpd=T,cex=2,inset=c(-0.1,0.02))
dev.off()
jpeg("top10_repeat_proportion.jpeg", width=1500,height=1200)
par(mai=c(3,3,1,3),mgp=c(13,0.9,0))
barplot(as.matrix(a),col=cols,main="Proportion of top 10 repeats",ylab="Order",xlab="Proportion of repeats",border=F,cex.axis=2,cex.names=1.5,cex.lab=1.5,horiz=T,las=1)
legend("topright",legend=rownames(a),fill=cols,bty="n",xpd=T,cex=2,inset=c(-0.16,0.02))
dev.off()
q()
n

