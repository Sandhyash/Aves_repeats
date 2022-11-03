## The input file is the repeat_summary_aves_ordered.txt file and the tree (aves.nwk) file
awk '{print $5}' repeat_summary_aves_ordered.txt |sed 's/[A-Z]\+/\n&/g'|sed 's/_$//g'|grep "[A-Z]"|sort|uniq -c|awk '{print $1,$2}' OFS="\t" > expansion_pic_species.txt
awk '{print $6}' repeat_summary_aves_ordered.txt |sed 's/[A-Z]\+/\n&/g'|sed 's/_$//g'|grep "[A-Z]"|sort|uniq -c|awk '{print $1,$2}' OFS="\t" > contraction_pic_species.txt
rm -f percentage_summary_pie.txt
for sp in `cat expansion_pic_species.txt contraction_pic_species.txt|cut -f2|sort -u`
do
echo "$sp"
exp=`awk -v s=$sp '$2==s{print $1}' expansion_pic_species.txt`
con=`awk -v s=$sp '$2==s{print $1}' contraction_pic_species.txt`
echo "$sp $exp $con"|awk '{print $1,$2,$3,$4=($2/($2+$3))*100,$5=($3/($2+$3))*100,$6=$2+$3}' OFS="\t" >> percentage_summary_pie.txt
done

## adding the pie charts on phylogeny
Rscript pie_chart.r
