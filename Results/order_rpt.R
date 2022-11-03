g=read.table("TROGONIFORMES_rpt_prop",header=F)
g1=g[order(g$V2,decreasing=T),]
write.table(g1,file="TROGONIFORMES_rpt_prop_ordered",quote=F,sep="\t",col.names=F,row.names=F)
