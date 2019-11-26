# This is used to remove overlapping peaks
:<<! First,the most significant peak is kept and any peak that directly overlaps with that significant peak is removed. Then, this process iterates to the next most significant peak and so on until all peaks have either been kept or removed due to direct overlap with a more significant peak. 
!

sample=K2
time=0
touch temporaryComplete.txt
split -l 10000 ${sample}_summits_blacklist_501.bed -d temporary_
ls -l | grep temporary_ | cut -d "_" -f 2 > temporaryIndex.txt
cat temporaryIndex.txt | while read line
do
{
	nrow=`wc -l temporary_$line | awk '{print $1}'`
	for ((i=1 ;i<$nrow; i=i+1))
	do
		Astart=`sed -n "${i}p" temporary_$line | awk '{print $2}'`
		Aend=`sed -n "${i}p" temporary_$line | awk '{print $3}'`
		Bstart=`sed -n "$[$i+1]p" temporary_$line | awk '{print $2}'`
		Bend=`sed -n "$[$i+1]p" temporary_$line | awk '{print $3}'`
		while [ `Rscript isoverlap.R ${Astart} ${Bstart} ${Aend} ${Bend}|awk '{print $2}'` == TRUE ]
		do
			Apvalue=`sed -n "${i}p" temporary_$line |awk '{print $5}'`
			Bpvalue=`sed -n "$[$i+1]p" temporary_$line | awk '{print $5}'`
			if [ `Rscript compare.R $Apvalue $Bpvalue |awk '{print $2}'` == TRUE ]
			then
				sed -i "$[$i+1]d" temporary_$line
			else
				sed -i "$[$i]d" temporary_$line
			fi
			Astart=`sed -n "${i}p" temporary_$line | awk '{print $2}'`
                	Aend=`sed -n "${i}p" temporary_$line | awk '{print $3}'`
                	Bstart=`sed -n "$[$i+1]p" temporary_$line | awk '{print $2}'`
                	Bend=`sed -n "$[$i+1]p" temporary_$line | awk '{print $3}'`
		done
		nrow2=` wc -l temporary_$line | awk '{print $1}'`
		if [ $i -ge $nrow2 ]
		then
	                echo "Complete" >> temporaryComplete.txt 
        	        break
      		fi
	done
} &
done
while [ `wc -l temporaryComplete.txt| awk '{print $1}'` != `wc -l temporaryIndex.txt| awk '{print $1}' ` ]
do
	
	echo "running"
	time=$[$time+30]
	sleep 30
done
cat temporary_* > ${sample}_summits_blacklist_501_rmoverlap.bed
rm temporary*
echo "Complete, time used:${time}s" 

# echo $sample
#  nrow=`wc -l ${sample}_summits_blacklist_501.bed | awk '{print $1}'` 
#  for ((i=1 ;i<$nrow; i=i+1))
#  do
#	Astart=`sed -n "${i}p" ${sample}_summits_blacklist_501.bed | awk '{print $2}'`
#	Aend=`sed -n "${i}p" ${sample}_summits_blacklist_501.bed | awk '{print $3}'`
#	Bstart=`sed -n "$[$i+1]p" ${sample}_summits_blacklist_501.bed | awk '{print $2}'`
#	Bend=`sed -n "$[$i+1]p" ${sample}_summits_blacklist_501.bed | awk '{print $3}'`

#	while [ `Rscript isoverlap.R ${Astart} ${Bstart} ${Aend} ${Bend}|awk '{print $2}'` == TRUE ]
#	do
# 		echo "do while"
#		Apvalue=`sed -n "${i}p" ${sample}_summits_blacklist_501.bed|awk '{print $5}'`
#		Bpvalue=`sed -n "$[$i+1]p" ${sample}_summits_blacklist_501.bed | awk '{print $5}'`
#		if [ `Rscript compare.R $Apvalue $Bpvalue |awk '{print $2}'` == TRUE ]
#		then
#			sed -i "$[$i+1]d" ${sample}_summits_blacklist_501.bed
#		#	echo "do then $i"
#		else
#			sed -i "$[$i]d" ${sample}_summits_blacklist_501.bed
#		#	echo "do else $[$i+1]"	
#		fi
#		Astart=`sed -n "${i}p" ${sample}_summits_blacklist_501.bed | awk '{print $2}'`
#	        Aend=`sed -n "${i}p" ${sample}_summits_blacklist_501.bed | awk '{print $3}'`
 #       	Bstart=`sed -n "$[$i+1]p" ${sample}_summits_blacklist_501.bed | awk '{print $2}'`
#        	Bend=`sed -n "$[$i+1]p" ${sample}_summits_blacklist_501.bed | awk '{print $3}'`
#	done
#	nrow2=` wc -l ${sample}_summits_blacklist_501.bed | awk '{print $1}'`
#	if [ $i -ge $nrow2 ]
#	then
#		echo "Break" 
#		break
#	fi
# done

































	 
