install.packages("devtools")    #必须用这个装wordcloud2，不然画不出来
devtools::install_github("lchiffon/wordcloud2",type="source") 
library(wordcloud2)
library(devtools)
install_github("lchiffon/wordcloud2")
figPath = system.file("examples/t.png",package = "wordcloud2")
wordcloud2(demoFreq, figPath = figPath, size = 1.5,color = "skyblue")
wordcloud2(demoFreq,figPath='C:/Users/Administrator/Desktop/123.png')
letterCloud(demoFreq, word = "R", size = 2)
library(tmcn)
cellword<-read.table("cellword.txt",sep="\t")
word_freq <- getWordFreq(string=unlist(cellword))
demoFreq
library(wordcloud2)  


