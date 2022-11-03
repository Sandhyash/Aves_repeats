library(phytools)
t1=read.tree("aves.nwk")
a=read.table("percentage_summary_pie.txt")
a$V4=a$V2/(a$V2+a$V3)
t=keep.tip(t1,a$V1)
b=data.frame()
for(sp in t$tip.label){
b=rbind(b,a[a$V1==sp,])}
rownames(b)=b$V1
b$V5=paste("(",b$V2,",",b$V3,")",sep="")
jpeg("aves_expan_con_pie.jpeg",width=20,height=20,units="in",res=300)
plot.phylo(t,use.edge.length=F,node.depth=2,label.offset=2.2,x.lim=30,cex=1.5,edge.width=2.5,main="Frequency of expanded and contracted genes in Aves clade")
tiplabels(pie=b$V4,piecol=c("orange","darkgreen"),cex=0.3)
tiplabels(b$V5,offset=1.2,cex=1.2,frame="none")
legend(x=23,y=3, legend = c("Genes with expanded repeats", "Genes with contracted repeats"), col=c("orange" , "darkgreen"),pch = 15, pt.cex = 3, cex = 1.5,  horiz = F,box.lwd = 0,box.col = "white",bg = "white")
dev.off()
q()
n


