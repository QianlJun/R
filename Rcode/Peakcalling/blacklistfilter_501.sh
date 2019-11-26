sample=K1
# filter peak summits by the ENCODE hg38 blacklist
cat ${sample}_summits.bed | awk '{print $1"\t"$3"\t"$3"\t""T""\t""T""\t"$4"\t"$5}' > ${sample}_summits.avinput
perl $annovar/table_annovar.pl ${sample}_summits.avinput $annovar/humandb/ -buildver hg38 -out $sample -remove -protocol ENCFF419RSJ -operation r -nastring .
sed -i '1d' ${sample}.hg38_multianno.txt

# after filtered, the peak summits were then extended by 250 bp on either side to final width of 501 bp 
paste ${sample}_summits.bed ${sample}.hg38_multianno.txt | awk '{if ($11==".") print }'| awk '{print $1"\t"$2-250"\t"$3+250"\t"$4"\t"$5}' > ${sample}_summits_blacklist_501.bed
