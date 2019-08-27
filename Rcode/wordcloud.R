setwd("D:/wanglab/workdata")
text<-readLines("cellword.txt" ) ##读取d取文档
text<-paste0(text,collapse = '') ##粘合每一行
text <- gsub("[^a-zA-Z]"," ",text) #除去非字母字符
text<- toupper(text)   #转化为大写
text<- strsplit(x= text,split =' ') #以空格为拆分符号
wordsFreq <-data.frame(table(text)) #建立data.frame
wordsFreq[,1]<-as.character(wordsFreq[,1]) #转化为character类型
wordsFreq<-wordsFreq[-which(nchar(wordsFreq[,1])<3),]  #过滤词频                            

wordsFreq<-wordsFreq[order(wordsFreq[,2],decreasing = T),] #给词频排序
write.table(wordsFreq,"wordsFreq.txt",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = ",") #保存词频文件

wordcloud<-read.table("wordsFreq.csv",sep = ",",header = TRUE)
library(wordcloud2)
figPath = system.file("examples/t.png",package = "wordcloud2")
wordcloud2(wordcloud, figPath = figPath, size = 1.5,color = "skyblue")
wordcloud2(wordcloud,figPath='C:/Users/Administrator/Desktop/333.png',color = "black",size=0.7)
letterCloud(demoFreq, word = "R", size = 2)

wordcloud<-read.table("wordcloudtest.csv",sep = ",",header = TRUE)
wordcloud2(wordcloud, size = 0.5,shape = 'star')  
