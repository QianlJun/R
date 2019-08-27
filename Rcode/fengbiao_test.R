# This the first part of the quality control process using data from fengbiao
raw<-read.table("G:/D/fengbiao_test/raw_countsHSC_sorted_revision1.txt",header=TRUE)
names(raw)
# Total number of reads>50000
cells.filtered<-c()
for (i in 2:ncol(raw)){
  if (sum(raw[i]) <= 50000){
    cells.filtered<-c(cells.filtered,i)
  }
}
cells.filtered
raw1<-raw[,-cells.filtered]
ncol(raw1)
# Number of genes detected (at least one mapped read) > 1000
names(raw1)
cells.filtered2<-c()
for (i in 2:ncol(raw1)){
  genes.number<-sum(raw1[,i]!=0)
  if (genes.number <= 1000){
    cells.filtered2<-c(cells.filtered2,i)
  }
}
cells.filtered2
raw2<-raw1[,-cells.filtered2]
ncol(raw2)
# Percentage of reads mapping to ERCC < 30%
ERCC.number<-length(grep(pattern = "ERCC",raw2$Gene_ID,value = T))
cells.filtered3<-c()
for (i in 2:ncol(raw2)){
  if (sum(raw2[1:ERCC.number,i])/sum(raw2[,i])>=0.3){
    cells.filtered3<-c(cells.filtered3,i)
  }
}
cells.filtered3
raw3<-raw2[,-cells.filtered3]
ncol(raw3)
# compute the number of cells passed all criteria
ncol(raw3)-1
# compute the number of oHSCs 
length(grep(pattern = "oHSCs",names(raw3),value = T))
# compute the number of dHSCs
length(grep(pattern = "dHSCs",names(raw3),value = T))
# output the data
write.csv(raw3,file="raw_countsHSC_sorted_revision_filitered.csv",quote = FALSE,row.names = FALSE)



# THis part is to find out highly variable genes beyond technical noise
library(DESeq)
library(genefilter)
library(statmod) 
raw_filtered<-read.table("G:/D/fengbiao_test/raw_countsHSC_sorted_revision_filitered.txt",header = TRUE,row.names = 1)
raw_filtered[1:10,1:2]
# Subset the table according to cell names:oHSCs, dHSCs
oH_countsAll<-raw_filtered[,grep(pattern="oHSCs",colnames(raw_filtered))]
dH_countsAll<-raw_filtered[,grep(pattern="dHSCs",colnames(raw_filtered))]

# Analysis of oHSCs 
# split the table into two sub-table, one with the ERCC spikes(countsERCC)
# one with the ENSM- genes(countsENSM)
geneTypes<-factor(c(EN="ENSM",ER="ERCC")[substr(rownames(raw_filtered),1,2)])
oH_countsENSM<-oH_countsAll[which(geneTypes=="ENSM"),]
oH_countsERCC<-oH_countsAll[which(geneTypes=="ERCC"),]
# Calculate size factors
oH_sfENSM<-estimateSizeFactorsForMatrix(oH_countsENSM)
oH_sfERCC<-estimateSizeFactorsForMatrix(oH_countsERCC)
rbind(oH_sfENSM,oH_sfERCC)
# Size-factor normalization
oH_nCountsENSM<-t(t(oH_countsENSM)/oH_sfENSM)
oH_nCountsERCC<-t(t(oH_countsERCC)/oH_sfERCC)
# Calculate the sample moments
oH_meansENSM<-rowMeans(oH_nCountsENSM)
oH_varsENSM<-rowVars(oH_nCountsENSM)
oH_cv2ENSM<-oH_varsENSM/oH_meansENSM^2

oH_meansERCC<-rowMeans(oH_nCountsERCC)
oH_varsERCC<-rowVars(oH_nCountsERCC)
oH_cv2ERCC<-oH_varsERCC/oH_meansERCC^2
# Fit technical noise
oH_minMeanForFit<-unname(quantile(oH_meansERCC[which(oH_cv2ERCC>0.3)],0.8))
oH_useForFit<-oH_meansERCC>=oH_minMeanForFit
oH_minMeanForFit
table(oH_useForFit)
fit<-glmgam.fit(cbind(a0=1,a1tilde=1/oH_meansERCC[oH_useForFit]),oH_cv2ERCC[oH_useForFit])
# Test for high variance
minBiolDisp<-0.5^2
xi<-mean(1/oH_sfERCC)
m<-ncol(oH_countsENSM)
psia1theta<-mean(1/oH_sfERCC)+(coefficients(fit)["a1tilde"]-xi)*mean(oH_sfERCC/oH_sfENSM)
cv2th<-coefficients(fit)["a0"]+minBiolDisp+coefficients(fit)["a0"]*minBiolDisp
testDenom<-(oH_meansENSM*psia1theta+oH_meansENSM^2*cv2th)/(1+cv2th/m)
p<-1-pchisq(oH_varsENSM*(m-1)/testDenom,m-1)
padj<-p.adjust(p,"BH")
table(padj<0.1)
# Select high variable genes
selectgenes<-names(subset(padj,padj<0.1))
oH_varibalegenes<-oH_nCountsENSM[selectgenes,]
nrow(oH_varibalegenes)
oH_varibalegenes[1:10,1:2]
# output data
setwd("G:/D/fengbiao_test/")
write.csv(oH_varibalegenes,file = "oHSCs_filtered_normalized.csv",quote = FALSE)

# Plot of results
plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 100 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) )
axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10" ,"100"), las=2 )
abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# Plot the plant genes, use a different color if they are highly variable
points( oH_meansENSM, oH_cv2ENSM, pch=20, cex=.2,
        col = ifelse( padj < .1, "#C0007090", "#70500040" ) )
# Add the technical noise fit, as before
xg <- 10^seq( -2, 6, length.out=1000 )
a0=1
lines( xg, coefficients(fit)["a1tilde"] / xg + a0, col="#FF000080", lwd=3 )
# Add a curve showing the expectation for the chosen biological CV^2 thershold
lines( xg, psia1theta/xg + coefficients(fit)["a0"] + minBiolDisp,
       lty="dashed", col="#C0007090", lwd=3 )
# Add the normalised ERCC points
points( oH_meansERCC, oH_cv2ERCC, pch=20, cex=1, col="#0060B8A0" )


# Analysis of dHSCs 
# split the table into two sub-table, one with the ERCC spikes(countsERCC)
# one with the ENSM- genes(countsENSM)
geneTypes<-factor(c(EN="ENSM",ER="ERCC")[substr(rownames(raw_filtered),1,2)])
dH_countsENSM<-dH_countsAll[which(geneTypes=="ENSM"),]
dH_countsERCC<-dH_countsAll[which(geneTypes=="ERCC"),]
dH_countsERCC[1:10,1:2]
write.csv(dH_countsAll,file="dHSCs_countsAll.csv",quote = FALSE)
# Calculate size factors
dH_sfENSM<-estimateSizeFactorsForMatrix(dH_countsENSM)
dH_sfERCC<-estimateSizeFactorsForMatrix(dH_countsERCC[,1:8])
dH_sfERCC
dH_countsERCC[,8]
rbind(dH_sfENSM,dH_sfERCC)



# Size-factor normalization
dH_nCountsENSM<-t(t(dH_countsENSM)/dH_sfENSM)
dH_nCountsERCC<-t(t(dH_countsERCC)/dH_sfERCC)
dH_nCountsERCC
# Calculate the sample moments
dH_meansENSM<-rowMeans(dH_nCountsENSM)
dH_varsENSM<-rowVars(dH_nCountsENSM)
dH_cv2ENSM<-dH_varsENSM/dH_meansENSM^2

dH_meansERCC<-rowMeans(dH_nCountsERCC)
dH_varsERCC<-rowVars(dH_nCountsERCC)
dH_cv2ERCC<-dH_varsERCC/dH_meansERCC^2
# Fit technical noise
dH_minMeanForFit<-unname(quantile(dH_meansERCC[which(dH_cv2ERCC>0.3)],0.8))
dH_useForFit<-dH_meansERCC>=dH_minMeanForFit
dH_minMeanForFit
table(dH_useForFit)
fit<-glmgam.fit(cbind(a0=1,a1tilde=1/dH_meansERCC[dH_useForFit]),dH_cv2ERCC[dH_useForFit])
# Test for high variance
minBiolDisp<-0.5^2
xi<-mean(1/dH_sfERCC)
m<-ncol(dH_countsENSM)
psia1theta<-mean(1/dH_sfERCC)+(coefficients(fit)["a1tilde"]-xi)*mean(dH_sfERCC/dH_sfENSM)
cv2th<-coefficients(fit)["a0"]+minBiolDisp+coefficients(fit)["a0"]*minBiolDisp
testDenom<-(dH_meansENSM*psia1theta+dH_meansENSM^2*cv2th)/(1+cv2th/m)
p<-1-pchisq(dH_varsENSM*(m-1)/testDenom,m-1)
padj<-p.adjust(p,"BH")
table(padj<0.1)
# Select high variable genes
selectgenes<-names(subset(padj,padj<0.1))
dH_varibalegenes<-dH_nCountsENSM[selectgenes,]
nrow(dH_varibalegenes)
dH_varibalegenes[1:10,1:2]
# output data
setwd("G:/D/fengbiao_test/")
write.csv(dH_varibalegenes,file = "dHSCs_filtered_normalized.csv",quote = FALSE)

# Plot of results
plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 100 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) )
axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10" ,"100"), las=2 )
abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# Plot the plant genes, use a different color if they are highly variable
points( dH_meansENSM, dH_cv2ENSM, pch=20, cex=.2,
        col = ifelse( padj < .1, "#C0007090", "#70500040" ) )
# Add the technical noise fit, as before
xg <- 10^seq( -2, 6, length.out=1000 )
a0=1
lines( xg, coefficients(fit)["a1tilde"] / xg + a0, col="#FF000080", lwd=3 )
# Add a curve showing the expectation for the chosen biological CV^2 thershold
lines( xg, psia1theta/xg + coefficients(fit)["a0"] + minBiolDisp,
       lty="dashed", col="#C0007090", lwd=3 )
# Add the normalised ERCC points
points( dH_meansERCC, dH_cv2ERCC, pch=20, cex=1, col="#0060B8A0" )


















# sfraw<-estimateSizeFactorsForMatrix(raw_filtered[1:10000,])
# sfraw
# help("estimateSizeFactorsForMatrix")
# head(countsENSM[1:100,])







# q<-read.table("raw_countsHSC_sorted_revision_filitered.txt",header = TRUE)
# 
# write.csv
# write.tab
# setwd("G:/D/fengbiao_test/")
# 
# a
# getwd()
# length(cells.filtered3)+length(cells.filtered2)+length(cells.filtered)
# sum(raw2[1:92,2])
# raw[raw$Gene_ID=="ERCC**",]
# grep("ERCC\\",raw2)
# raw1[,79]
# sum(raw1[,84]!=0)
# raw1[,]
# names(raw1)
# raw1[,2!=0]
# raw1[,2]
# head(raw1[,2])
# sum(raw1[,2]==0)
# sum(raw1[,2]!=0)
# x="hello world"
# grep(pattern = "hello",x)
# grep(pattern = "l",x)
# grep(pattern = "h",x,value = T)
# a<-grep(pattern = "ERCC",raw2[,1],value = T)
# raw2[,1]
# length(a)
# sum(raw2[1:ERCC.number,79])/sum(raw2[,79])<0.3
# sum(raw2[1:ERCC.number,2])
# sum(raw2[,2])
# sum(raw2[1:ERCC.number,i])/sum(raw2[,i]>=0.3