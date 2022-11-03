library("UpSetR")
a<-read.table("orderwise_preabs_final",header=T)
b<-read.table("count_totalgenes_orderwise")
## adding the column for the proortions of genes in all order
allsets<-colnames(a[,-1])
allset<-data.frame(order=allsets)
for (i in allset$order){
allset[allset$order==i,2]<-sum(a[,(grep(i, colnames(a)))])/(b[b$V1==colnames(a[i]),2])}
colnames(allset)[2]<-c("Proportion_of_genes")
pdf("upset_plot_orderwise.pdf",width=14,height=7)
par(mar=c(0.1,12,0.1,0.1))
upset(a,sets=b$V1,order.by="freq", matrix.color=c("blue"),main.bar.color="darkgreen",mainbar.y.label="Number of genes",decreasing=T,keep.order=T,point.size=2,text.scale=c(1.3,1.3,1,0.8,0.6,0.75), set.metadata = list(data = allset, plots = list(list(type = "hist",column = "Proportion_of_genes", assign = 20,colors="plum"))), sets.bar.color=c("saddlebrown"),nintersects = 24,sets.x.label="Number of genes")
dev.off()
jpeg("upset_plot_orderwise.jpeg",width=1500,height=1500)
par(mai=c(0.1,1,0.1,0.1))
upset(a,sets=b$V1,order.by="freq", matrix.color=c("blue"),main.bar.color="darkgreen",mainbar.y.label="Number of genes",decreasing=T,keep.order=T,point.size=3,text.scale=c(2.3,2.3,2,2,1.9,1.75), set.metadata = list(data = allset, plots = list(list(type = "hist",column = "Proportion_of_genes", assign = 15,colors="plum",text.scale=3))), sets.bar.color=c("saddlebrown"),nintersects = 24,sets.x.label="Number of genes")
dev.off()
q()
n

