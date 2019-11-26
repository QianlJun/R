# use to combine reproducible peaks in all sample
# sort according to chromosome and position

cat *_normalized_reproducible.bed > reproducible_all.bed
sort -k1,1 -k2,2n reproducible_all.bed > reproducible_all_sort.bed

