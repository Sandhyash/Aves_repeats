##This script makes the stacked barplot, upset plot, and simpson diversity plot
#All the major analyses will be performed on aves_morethanseventy_combined
#species_orders.txt contains names of species and orders
#birds_orf_list file contains list of all the genes appearing in at least one bird species #15042

#########Making the figures###################
###1. Stacked bar-plot
awk '($8-$7)>3 && ($9/($8-$7))>0.7 && length($11)<5{print $3"\t"$2"\t"$11}' aves_morethanseventy_combined|sed 's/\([^[:blank:]]\)\([[:upper:]]\)/\1\ \2/2'|awk '{print $2"\t"$3"\t"$4}' > species_rpts

while read spt
do
echo $spt
ord=`echo $spt|awk '{print $1}'`
sp=`echo $spt|awk '{print $2}'`
grep "\b$sp\b" species_rpts|awk -v o=$ord '{print o"\t"$1"\t"$3}' >> "$ord"_rpt
done < species_orders.txt

for ord in `awk '{print $1}' species_orders.txt|sort -u`
do
echo $ord
len=`wc -l "$ord"_rpt|awk '{print $1}'`
awk '{print $3}' "$ord"_rpt|sort|uniq -c|awk -v l=$len '{print $2"\t"$1/l}'|sort -k2nr > "$ord"_rpt_prop
done
###The above sort command does not sort the numerical values if they are written in scientific format; We have to use R#########
for ord in `awk '{print $1}' species_orders.txt|sort -u`
do
echo $ord
echo 'g=read.table("'"$ord"'_rpt_prop",header=F)' > order_rpt.R
echo 'g1=g[order(g$V2,decreasing=T),]' >> order_rpt.R
echo 'write.table(g1,file="'"$ord"'_rpt_prop_ordered",quote=F,sep="\t",col.names=F,row.names=F)' >> order_rpt.R
Rscript order_rpt.R
done
######top 20 has a lot of amino acids; now doing it with top ten###
head -n10 *_rpt_prop_ordered|grep -v "="|awk '{print $1}'|sort -u|sed '/^$/d' > all_top10_rpts
echo -e "Repeats\tPSITTACIFORMES\tPASSERIFORMES\tACCIPITRIFORMES\tPELECANIFORMES\tSPHENISCIFORMES\tGALLIFORMES\tANSERIFORMES\tTINAMIFORMES" > top10_rpt_prop_combined
while read rpt
do
echo $rpt
echo "$rpt" > temp.txt
for ord in PSITTACIFORMES PASSERIFORMES ACCIPITRIFORMES PELECANIFORMES SPHENISCIFORMES GALLIFORMES ANSERIFORMES TINAMIFORMES
do
awk -v r=$rpt '$1==r{print $2}' "$ord"_rpt_prop_ordered >> temp.txt
done
sed -z 's/\n/\t/g' temp.txt|sed 's/\(.*\)\t/\1\n/g' > temp2.txt
cat top10_rpt_prop_combined temp2.txt >> top10_rpt_prop_combined
done < all_top10_rpts
rm temp*
##check the number of unique repeats in your top10 file: all_top10_rpts
#That will decide the position of "Remaining" in the below R script
##In R##
#R
#a<-read.table("top10_rpt_prop_combined",header=T)
#a$Repeats<-as.character(a$Repeats)
#a[14,1]<-"Remaining" #add this row at the end
#for(i in 2:length(a)){
#a[14,i]<-1-sum(na.omit(a[,i]))
#}
#a<-data.frame(a,row.names=1)
#cols=c("yellowgreen","violetred","violet","turquoise1","tomato","thistle","tan","steelblue1","springgreen3","slateblue2","salmon","purple","plum1","saddlebrown")
#pdf("top10_repeat_proportion.pdf",width=30,height=20)
#par(mai=c(3,3,1,3),mgp=c(13,0.9,0))
#barplot(as.matrix(a),col=cols,main="Proportion of top 10 repeats",ylab="Order",xlab="Proportion of repeats",border=F,cex.axis=2,cex.names=1.5,cex.lab=1.5,horiz=T,las=1)
#legend("topright",legend=rownames(a),fill=cols,bty="n",xpd=T,cex=2,inset=c(-0.1,0.02))
#dev.off()
#jpeg("top10_repeat_proportion.jpeg", width=1500,height=1200)
#par(mai=c(3,3,1,3),mgp=c(13,0.9,0))
#barplot(as.matrix(a),col=cols,main="Proportion of top 10 repeats",ylab="Order",xlab="Proportion of repeats",border=F,cex.axis=2,cex.names=1.5,cex.lab=1.5,horiz=T,las=1)
#legend("topright",legend=rownames(a),fill=cols,bty="n",xpd=T,cex=2,inset=c(-0.16,0.02))
#dev.off()
#q()
#n
Rscript stacked.r

echo "stacked-barplot made"
################################################Stacked barplot ends here####################################
###2. Simpson's diversity index
for rptfile in `ls *_rpt`
do
ord=`echo $rptfile|sed 's/_rpt//g'`
echo $ord
awk '{print $3}' $rptfile|sort|uniq -c|awk '{print $2,$1}' OFS="\t" > "$ord"_rpt_count
done
rm -f simp_d_orderwise.txt

for cfiles in `ls *_rpt_count`
do
ord=`echo $cfiles|sed 's/_rpt_count//g'`
echo 'a=read.table("'"$ord"'_rpt_count",header=F)' > simp.r
echo 'a$V3=(a$V2*(a$V2 -1))/(sum(a$V2)*(sum(a$V2)-1))' >> simp.r
echo 'd=1-sum(a$V3)' >> simp.r
echo 'e=paste("'$ord'",d,sep="\t")' >> simp.r
echo 'write.table(e,file="simp_d_orderwise.txt",append=T,quote=F,sep="\t",col.names=F,row.names=F)' >> simp.r
Rscript simp.r
done
##########Plotting barplot of Simpson's diversity index########
##In R##
Rscript simpson_diversity.r
#a=read.table("simp_d_orderwise.txt",header=F)
#library(stringr)
#a$V1=str_to_title(a$V1)
#jpeg("simpson_diversity_aves_repeats.jpeg",width=23.2,height=13,units="in",res=300)
#par(mai=c(2.6,2.6,0.5,0.2),mgp=c(11,0.2,0))
#barplot(a$V2,names=a$V1,ylab="",ylim=c(0,1),xlab="",col="#F7B37E",main="Simpson diversity of repeats across different orders",las=3,cex.axis=1,cex.lab=1.5)
#title(xlab="Orders",line=8,cex.lab=2)
#title(ylab="Simpson diversity",cex.lab=2,line=2)
#dev.off()
echo "simpsons deversity plot made"
#################################################Simpson's diversity ends here######################################
###3. UpSet plot
awk '($8-$7)>3 && ($9/($8-$7))>0.7 && length($11)<5{print $3"\t"$2"\t"$11}' aves_morethanseventy_combined|sed 's/\([^[:blank:]]\)\([[:upper:]]\)/\1\ \2/2'|cut -f2 -d" " > sp_gene_rpt
while read j
do
echo $j
sp=`echo $j|awk '{print $2}'`
ord=`echo $j|awk '{print $1}'`
grep "\b$sp\b" sp_gene_rpt|awk -v o=$ord '{print o,$0}' OFS="\t" >> ord_sp_gene_rpt
done < species_orders.txt

####making a list of all the repeat-containing genes present in the ord_sp_gene_rpt file
awk '{print $3}' ord_sp_gene_rpt |sort -u > rpt_genelist
for ord in `awk '{print $1}' ord_sp_gene_rpt |sort -u|sed -z 's/\n/ /g'`
do
grep "$ord" ord_sp_gene_rpt > "$ord"_generp.txt
done
rm -f *_preabs

###We are making the prsence-absence matrix in a orderwise manner for each repeat-containing gene
for ord in `awk '{print $1}' ord_sp_gene_rpt |sort -u|sed -z 's/\n/ /g'`
do
echo $ord
while read gene
do
num=`grep -c "\b$gene\b" "$ord"_generp.txt`
if [ $num -ge 1 ]
then
echo "1" >> "$ord"_preabs
else
echo "0" >> "$ord"_preabs
fi
done < rpt_genelist
echo "$ord done"
done

paste rpt_genelist ACCIPITRIFORMES_preabs -d "\t" > orderwise_preabs
for ord in `awk '{print $1}' ord_sp_gene_rpt |sort -u|sed '1d'|sed -z 's/\n/ /g'`
do
paste orderwise_preabs "$ord"_preabs -d "\t" > temp
mv temp orderwise_preabs
done
paste <(awk '{print $1}' ord_sp_gene_rpt |sort -u|sed '1 i\gene'|sed -z 's/\n/ /g') orderwise_preabs -d "\n"|sed '/^$/d' > orderwise_preabs_final

###Now making the proportion file of the repeat containing genes
##firstly, total number of genes in Aves
##we have a file birds_orf_list which has names of all the genes appearing in at least one bird species
##15042 genes
##Now for each order we have to check the proportion of genes containing a repeat
##we will make use of the list containing species names and order name (species_orders.txt)

#rm -f gene_species.txt
#for file in `ls all_protein_orfs/*birds_list_orf.fa`
#do
#echo $file
#gene=`echo $file|sed 's#all_protein_orfs/##g'|cut -d"_" -f1`
#grep ">" $file|sed -e 's/>//g' -e 's/\([^[:blank:]]\)\([[:upper:]]\)/\1\ \2/2'|cut -f2 -d" " |sed 's/\n/ /g'|awk -v g=$gene '{print g,$0}' >> gene_species.txt
#done
###sorting the above made file according to the orders
while read j
do
echo $j
ord=`echo $j|awk '{print $1}'`
sp=`echo $j|awk '{print $2}'`
grep "\b$sp\b" gene_species.txt >> "$ord"_protein_list
done < species_orders.txt

for file in `ls *_protein_list`
do
ord=`echo $file|sed 's/_protein_list//g'`
awk '{print $1}' $file|sort -u > "$ord"_unique_gene
done
##The file *_unique_gene contains unique protein coding genes in each order

###making the proportion repeats file for each order
###We will include the proportion in the UpSet plot
##The proportion is defined as the number of genes with repeats to the number of protein coding genes in that order
rm -f proportion_genes_rpt_orderwise
for file in `ls *_unique_gene`
do
echo $file
totgene=`wc -l $file|awk '{print $1}'`
ord=`echo $file|sed 's/_unique_gene//g'`
rpfile=`ls "$ord"_generp.txt`
rpgenenum=`awk '{print $3}' $rpfile|sort -u|wc -l|awk '{print $1}'`
echo $rpgenenum $totgene|awk -v o=$ord '{print o,$1/$2}' OFS="\t" >> proportion_genes_rpt_orderwise
done
############making the above file in the count also
rm -f count_totalgenes_orderwise
for file in `ls *_unique_gene`
do
echo $file
totgene=`wc -l $file|awk '{print $1}'`
ord=`echo $file|sed 's/_unique_gene//g'`
echo $totgene|awk -v o=$ord '{print o,$1}' OFS="\t" >> count_totalgenes_orderwise
done

###making number of repeat-containing genes in each order
rm -f count_genes_rpt_orderwise
for file in `ls *_unique_gene`
do
echo $file
ord=`echo $file|sed 's/_unique_gene//g'`
rpfile=`ls "$ord"_generp.txt`
rpgenenum=`awk '{print $3}' $rpfile|sort -u|wc -l|awk '{print $1}'`
echo $rpgenenum|awk -v o=$ord '{print o,$1}' OFS="\t" >> count_genes_rpt_orderwise
done
#################in R##########
#library("UpSetR")
#a<-read.table("orderwise_preabs_final",header=T)
#b<-read.table("count_totalgenes_orderwise")
### adding the column for the proortions of genes in all order
#allsets<-colnames(a[,-1])
#allset<-data.frame(order=allsets)
#for (i in allset$order){
#allset[allset$order==i,2]<-sum(a[,(grep(i, colnames(a)))])/(b[b$V1==colnames(a[i]),2])}
#colnames(allset)[2]<-c("Proportion_of_genes")
#pdf("upset_plot_orderwise.pdf",width=14,height=7)
#par(mar=c(0.1,12,0.1,0.1))
#upset(a,sets=b$V1,order.by="freq", matrix.color=c("blue"),main.bar.color="darkgreen",mainbar.y.label="Number of genes",decreasing=T,keep.order=T,point.size=2,text.scale=c(1.3,1.3,1,0.8,0.6,0.75), set.metadata = list(data = allset, plots = list(list(type = "hist",column = "Proportion_of_genes", assign = 20,colors="plum"))), sets.bar.color=c("saddlebrown"),nintersects = 24,sets.x.label="Number of genes")
#dev.off()
#jpeg("upset_plot_orderwise.jpeg",width=1500,height=1500)
#par(mai=c(0.1,1,0.1,0.1))
#upset(a,sets=b$V1,order.by="freq", matrix.color=c("blue"),main.bar.color="darkgreen",mainbar.y.label="Number of genes",decreasing=T,keep.order=T,point.size=3,text.scale=c(2.3,2.3,2,2,1.9,1.75), set.metadata = list(data = allset, plots = list(list(type = "hist",column = "Proportion_of_genes", assign = 15,colors="plum",text.scale=3))), sets.bar.color=c("saddlebrown"),nintersects = 24,sets.x.label="Number of genes")
#dev.off()
#q()
#n
Rscript upset.r
echo "upset plot made"
################################################################UpSet plot ends here######################################################


