# Using MACS2 to callpeak
sample=K1
macs2 callpeak -t ${sample}_sort_dedup.bam -f BAM -g hs --outdir ~/ATAC-seq/data/ANCBJ170581_PM-BJ170581-04_AHWHHNCCXY_2018-12-10/Rawdata/aligned_deduplicated_q10/$sample -n $sample --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01
