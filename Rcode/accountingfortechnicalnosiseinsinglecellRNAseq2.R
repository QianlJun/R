# This part is to find out highly variable genes beyond technical noise
library(DESeq)
library(genefilter)
library(statmod) 
library( org.At.tair.db )

raw_filtered<-read.csv("raw_countsHSC_sorted_revision_filitered.csv",header = TRUE,row.names = 1)
raw_filtered[1:10,1:2]
# Filter genes without any read map to them 
fullCountTable<-raw_filtered[rowSums(raw_filtered)!=0,]
fullCountTable[1:10,1:2]
# fullCountTable<-raw_filtered
# log(x+1)-transform
# fullCountTable_log<-log(fullCountTable+0.1)
fullCountTable_log<-(fullCountTable+0.1)
fullCountTable_log[1:10,1:2]
# split the table into two sub-table,one with the ERCC spikes,one with the mouse genes
geneTypes<-factor(c(EN="Mmus",ER="ERCC")[substr(rownames(fullCountTable_log),1,2)])
countsMmus<-fullCountTable_log[which(geneTypes=="Mmus"),]
countsERCC<-fullCountTable_log[which(geneTypes=="ERCC"),]
countsMmus[1:10,1:2]
countsERCC[1:10,1:2]


sfERCC<-estimateSizeFactorsForMatrix(countsERCC)
sfMmus<-estimateSizeFactorsForMatrix(countsMmus)
rbind(sfMmus,sfERCC)
nCountsMmus<-t(t(countsMmus)/sfMmus)
nCountsERCC<-t(t(countsERCC)/sfERCC)

meansERCC <- rowMeans( nCountsERCC )
varsERCC <- rowVars( nCountsERCC )
cv2ERCC <- varsERCC / meansERCC^2
minMeanForFit<-unname(quantile(meansERCC[which(cv2ERCC>.3)],.8))
minMeanForFit

useForFit<-meansERCC>=minMeanForFit
fit<-glmgam.fit(cbind(a0=1,a1tilde=1/meansERCC[useForFit]),cv2ERCC[useForFit])
fit$coefficients

xi<-mean(1/sfERCC)
a0<-unname(fit$coefficients["a0"])
a1<-unname(fit$coefficients["a1tilde"]-xi)
c(a0,a1)


meansMmus <- rowMeans( nCountsMmus )
varsMmus <- rowVars( nCountsMmus )
cv2Mmus <- varsMmus / meansMmus^2
psia1theta<-mean(1/sfMmus)+a1*mean(sfERCC/sfMmus)
minBiolDisp<-0.5^2

m<-ncol(countsMmus)
cv2th<-a0+minBiolDisp+a0*minBiolDisp
testDenom<-(meansMmus*psia1theta+meansMmus^2*cv2th)/(1+cv2th/m)
p<-1-pchisq(varsMmus*(m-1)/testDenom,m-1)
padj<-p.adjust(p,"BH")
sig<-padj<.1
sig[is.na(sig)]<-FALSE
table(sig)

log2RelExprMmus<-log2(nCountsMmus/meansMmus)
highVarTable<-data.frame(row.names=NULL,geneID=rownames(countsMmus)[sig],
                         meanNormCount=meansMmus[sig],
                         strongest=factor(colnames(log2RelExprMmus)[
                           apply(log2RelExprMmus[sig,],1,which.max)]),
                         log2RelExprMmus[sig,],
                         check.names=FALSE)
highVarTable[1:10,1:5]
write.csv(highVarTable,file="highly_variant_genes_Mmus.csv",row.names = FALSE)
