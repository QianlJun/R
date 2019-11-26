# normalize peak significance scores and select reproducible peak(score per million value >=5)
# method: each(pvalue)/(sum(pvalue)/1000000)
args=commandArgs(T)
data<-read.table(args[2])
#head(data)
data$V6<-data$V5/(sum(data$V5)/1000000)
data_normalized5<-data[data$V6>=5,]
data_normalized_reproducible<-data_normalized5[,c(1,2,3,4,6)]
filename_nor<-paste(args[1],"_normalized.bed",sep="")
filename_rep<-paste(args[1],"_normalized_reproducible.bed",sep="")
write.table(data,file=filename_nor,sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(data_normalized_reproducible,file=filename_rep,sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)



