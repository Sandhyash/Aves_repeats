####This script will plot everything
args = commandArgs(trailingOnly=TRUE)
if (length(args)<6) {
  stop("Please provide arguments like Rscript repeatplotterwitheverything.r inputfile treefile cladename genename windowsize codemlresult ", call.=FALSE)
}
library(ape)
library(stringr)
library(TreeTools)
library(tidytree)
a1=read.table(args[1])
t1=read.tree(args[2])
colnames(a1)=c("clade","gene","species","rpt","unalnstart","unalnend","alnstart","alnend","alnlength")
clade=tolower(args[3])
cname=str_to_title(clade)
gene=toupper(args[4])
proximval=as.numeric(args[5])
ssfile=read.table(args[6])
colnames(ssfile)=c("clade","gene","aasite","aa","pval")
if(proximval>20){proximval=20}
colpal=data.frame(aa = c("A","B","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","U","V","W","X","Y","Z"), col=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","orange","cyan","violet"))
if( any(a1$gene==gene) && dim(a1[a1$gene==gene,])[1] > 3 && any(a1$clade==clade) && dim(a1[a1$clade==clade,])[1] > 3 |any(a1$gene==gene) && dim(a1[a1$gene==gene,])[1] > 3 && cname=="All") {
if(cname=="All"){a=a1
ss=ssfile[ssfile$gene==gene,] } else {
a<-a1[a1$clade==clade & a1$gene==gene,]
ss=ssfile[ssfile$clade==clade & ssfile$gene==gene,]}
final=a
t<-keep.tip(t1,unique(c(a$species)))
pt=Preorder(t)
heightparam=max(((length(t$tip.label))*0.435)+1.70,14)
species=as.data.frame(t$tip.label)
colnames(species)=c("spnames")
final$rplen=(final$unalnend - final$unalnstart ) + 1
aas=data.frame(unique(as.data.frame(strsplit(paste(final$rpt,collapse=""),split="")))[,1])
colnames(aas)=c("aa")
newcols=merge(aas,colpal,by="aa")
jpeg(paste('Plot_for_',gene,'_in_',cname,'withpic.jpeg',sep=''),width=28,height=heightparam,units="in",res=300)
par(mai=c(0.4,0,0.4,0),oma=c(0,0.4,0,0.2),xpd=T)
layout(matrix(c(1,2,2,2),nrow=1,byrow=T))
###Firstly, plotting the species tree for the gene
phylolim=max(((length(t$tip.label))/5) + 6,14)
plot.phylo(t,use.edge.length=F,node.depth=2,label.offset=0.1,x.lim=phylolim,cex=2,edge.width=4.5)
ymax=2*length(t$tip.label)
plot(1, type="n", xlab="", ylab="", xlim=c(0, max(final$alnlength)), ylim=c(0, ymax-2),axes=F)
###Plotting the recatangular boxes on the right side of the plot
for(j in seq(0,length(t$tip.label)-1,1)){
y1=(2*j)-0.5
y2=(2*j)+0.5
rect(0,y1,max(final$alnlength),y2)}
###Plotting codeml results in the rectangles
if(dim(ss)[1] > 0){
for(l in unique(ss$aasite)){
x1=3*l
y1=-0.5
x2=(3*l)+2
y2=ymax-1.5
rect(x1,y1,x2,y2,border="transparent",col=rgb(0.5,0,0,alpha=0.3)) } }
###Plotting the repeats at their specific coordinates
for (z in 1:dim(final)[1]){
b1=final[z,]
spnaam=b1[1,"species"]
###get the number at which the species is written
j=which(species$spnames==spnaam)
x1=b1$alnstart
x2=b1$alnstart + b1$rplen
s=as.character(b1$rpt)
#print(s)
y1=((2*j)-0.5)-2
y2=((2*j)+0.5)-2
if(nchar(s)>4){
for(i in 1:nchar(s)){
if(i<=nchar(s)-1){
aanow=data.frame(strsplit(s,split=""))[i,]
meracol=newcols[newcols$aa==aanow,2]
rect(x1,y1,x2,y2,brder="transparent",col=meracol,lty="dashed",density=40,angle=45*i)}
else {aanow=data.frame(strsplit(s,split=""))[5,]
meracol=newcols[newcols$aa==aanow,2]
rect(x1,y1,x2,y2,border=meracol,lwd=3)
}
}
} else {for(i in 1:nchar(s)){
aanow=data.frame(strsplit(s,split=""))[i,]
meracol=newcols[newcols$aa==aanow,2]
rect(x1,y1,x2,y2,border="transparent",col=meracol,lty="dashed",density=40,angle=45*i)
}
}
rplabel=paste(b1$rpt,' (',b1$rplen,') ',sep='')
text((x1+x2)/2,y2+0.3,labels=rplabel)
}
ld=data.frame(unique(as.data.frame(strsplit(paste(final$rpt,collapse=""),split="")))[,1])
colnames(ld)=c("aa")
aminocol<-merge(ld,colpal,by="aa")
legend("bottomright",legend=aminocol$aa,col=aminocol$col,horiz=T,fill=aminocol$col,inset=c(0.042,-0.01),bty="n",cex=2.4)
genenm=final$gene[1]
ti=substitute('Repeat length distribution of'~italic(genenaam)~'gene across different species in'~cladename~'clade',list(genenaam=genenm,cladename=cname))
mtext(ti,side=3,cex=2,adj=0.5,outer=T,line=-4)
m=final
if (length(unique(a$alnend)) < 2) { en =as.numeric(unique(a$alnend))
x = m[sqrt((m$alnend-en)**2)<=proximval, ]
rrp=x[1,"rpt"]
if(dim(x)[1]>=4){ 
ps=x[c("species","rplen")]
trp=Preorder(keep.tip(t,as.character(ps$species)))
trpsp=trp$tip.label
psrp=data.frame()
for(sps in trpsp){
df1=ps[ps$species==sps,]
psrp=rbind(psrp,df1)
}
pspname=psrp
psrplen<-data.frame(psrp[,-1])
rownames(psrplen)=as.character(pspname$species)
colnames(psrplen)<-c("repeat_length")
pic.rplen<-pic(psrplen$repeat_length,trp)
bigval=quantile(pic.rplen,0.99)
smallval=quantile(pic.rplen,0.01)
allsigs=as.data.frame(pic.rplen[pic.rplen> bigval | pic.rplen< smallval])
tbltrp=as_tibble(trp)
for(fz in as.numeric(rownames(allsigs))){
ch=child(tbltrp,fz)
if (length(as.numeric(ch$node))==2){
ns1=as.numeric(ch$node)[1]
ns2=as.numeric(ch$node)[2]
spns1=Subtree(trp,ns1)$tip.label
spns2=Subtree(trp,ns2)$tip.label
if (length(spns1) > 2 && length(spns2) > 2){
dfns1=x[x$species %in% spns1,]
dfns2=x[x$species %in% spns2,]
if( mean(dfns1$rplen) > mean(dfns2$rplen) ) { colvals=c("red","blue")} else { colvals=c("blue","red")}
ns1x1=min(dfns1$alnstart) - 2
ns1x2=max(dfns1$alnstart + dfns1$rplen) + 2
ns1y1=100000
for (spnaam in dfns1$species){
ns1y1=min(ns1y1,which(species$spnames==spnaam)) }
ns1y2=0
for (spnaam in dfns1$species){
ns1y2=max(ns1y2,which(species$spnames==spnaam)) }
y1=(2*(ns1y1) - 1)-2
y2=(2*(ns1y2)+1)-2
rect(ns1x1,y1,ns1x2,y2,density=0,border=colvals[1],lwd=2.2)
ns2x1=min(dfns2$alnstart) - 2
ns2x2=max(dfns2$alnstart + dfns2$rplen) + 2
ns2y1=100000
for (spnaam in dfns2$species){
ns2y1=min(ns2y1,which(species$spnames==spnaam)) }
ns2y2=0
for (spnaam in dfns2$species){
ns2y2=max(ns2y2,which(species$spnames==spnaam)) }
y1=(2*(ns2y1) - 1)-2
y2=(2*(ns2y2)+1)-2
rect(ns2x1,y1,ns2x2,y2,density=0,border=colvals[2],lwd=2.2)
fortbl=paste(clade,gene,rrp,paste0(as.vector(dfns1$species),collapse='_'),paste0(as.vector(dfns2$species),collapse='_'),mean(dfns1$rplen),mean(dfns2$rplen),fz,collapse="\t")
tabl=paste(gene,"repeat_summary.txt",sep="")
write.table(fortbl,file=tabl,append=T,sep='\t', quote = F,row.names = F, col.names = F)
}
}
}
}
if(dim(m)[1]==0){break}
} else { nbreaks=length(seq(from=min(final$alnstart),to=max(final$alnend),by=proximval))
h=hist(final$alnend,breaks=nbreaks,plot=F)
bi=data.frame(h$mids,h$counts)
bis=bi[bi$h.counts!=0,]
bis=bis[order(-bis$h.counts),]
for( en in as.list(bis$h.mids)){
x = m[sqrt((m$alnend-en)**2)<=proximval, ]
rrp=x[1,"rpt"]
m=m[sqrt((m$alnend-en)**2)>proximval, ]
if(dim(x)[1]>=4){ 
ps=x[c("species","rplen")]
trp=Preorder(keep.tip(t,as.character(ps$species)))
trpsp=trp$tip.label
psrp=data.frame()
for(sps in trpsp){
df1=ps[ps$species==sps,]
psrp=rbind(psrp,df1)
}
pspname=psrp
psrplen<-data.frame(psrp[,-1])
rownames(psrplen)=as.character(pspname$species)
colnames(psrplen)<-c("repeat_length")
pic.rplen<-pic(psrplen$repeat_length,trp)
bigval=quantile(pic.rplen,0.99)
smallval=quantile(pic.rplen,0.01)
allsigs=as.data.frame(pic.rplen[pic.rplen> bigval | pic.rplen< smallval])
tbltrp=as_tibble(trp)
for(fz in as.numeric(rownames(allsigs))){
ch=child(tbltrp,fz)
if (length(as.numeric(ch$node))==2){
ns1=as.numeric(ch$node)[1]
ns2=as.numeric(ch$node)[2]
spns1=Subtree(trp,ns1)$tip.label
spns2=Subtree(trp,ns2)$tip.label
if (length(spns1) > 2 && length(spns2) > 2){
dfns1=x[x$species %in% spns1,]
dfns2=x[x$species %in% spns2,]
if( mean(dfns1$rplen) > mean(dfns2$rplen) ) { colvals=c("red","blue")} else { colvals=c("blue","red")}
ns1x1=min(dfns1$alnstart) - 2
ns1x2=max(dfns1$alnstart + dfns1$rplen) + 2
ns1y1=100000
for (spnaam in dfns1$species){
ns1y1=min(ns1y1,which(species$spnames==spnaam)) }
ns1y2=0
for (spnaam in dfns1$species){
ns1y2=max(ns1y2,which(species$spnames==spnaam)) }
y1=(2*(ns1y1) - 1)-2
y2=(2*(ns1y2)+1)-2
rect(ns1x1,y1,ns1x2,y2,density=0,border=colvals[1],lwd=2.2)
ns2x1=min(dfns2$alnstart) - 2
ns2x2=max(dfns2$alnstart + dfns2$rplen) + 2
ns2y1=100000
for (spnaam in dfns2$species){
ns2y1=min(ns2y1,which(species$spnames==spnaam)) }
ns2y2=0
for (spnaam in dfns2$species){
ns2y2=max(ns2y2,which(species$spnames==spnaam)) }
y1=(2*(ns2y1) - 1)-2
y2=(2*(ns2y2)+1)-2
rect(ns2x1,y1,ns2x2,y2,density=0,border=colvals[2],lwd=2.2)
fortbl=paste(clade,gene,rrp,paste0(as.vector(dfns1$species),collapse='_'),paste0(as.vector(dfns2$species),collapse='_'),mean(dfns1$rplen),mean(dfns2$rplen),fz,collapse="\t")
tabl=paste(gene,"repeat_summary.txt",sep="")
write.table(fortbl,file=tabl,append=T,sep='\t', quote = F,row.names = F, col.names = F)
}
}
if(dim(m)[1]==0){break}
}}}}
dev.off()
}else {
print("Either the clade/gene name is missing/wrong or number of species < 4 or gene is missing in gene_length file")
break }
