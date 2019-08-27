# This the first part of the quality control process using data from fengbiao
setwd("G:/D/fengbiao_test/")
raw<-read.table("raw_countsHSC_sorted_revision.txt",header=TRUE)
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



# This part is to find out highly variable genes beyond technical noise
library(DESeq)
library(genefilter)
library(statmod) 
raw_filtered<-read.csv("raw_countsHSC_sorted_revision_filitered.csv",header = TRUE,row.names = 1)
raw_filtered[1:10,1:2]
# Filter genes without any read map to them 
fullCountTable<-raw_filtered[rowSums(raw_filtered)!=0,]
fullCountTable[1:10,1:2]

# log(x+0.1)-transform(skip)
# fullCountTable_log<-log(fullCountTable+0.1)
fullCountTable_log<-fullCountTable+0.1
fullCountTable_log[1:10,1:2]
# split the table into two sub-table,one with the ERCC spikes,one with the mouse genes
geneTypes<-factor(c(EN="Mmus",ER="ERCC")[substr(rownames(fullCountTable_log),1,2)])
countsMmus<-fullCountTable_log[which(geneTypes=="Mmus"),]
countsERCC<-fullCountTable_log[which(geneTypes=="ERCC"),]
countsMmus[1:10,1:2]
countsERCC[1:10,1:2]

# normalized by size factors
sfMmus<-estimateSizeFactorsForMatrix(countsMmus)
sfERCC<-estimateSizeFactorsForMatrix(countsERCC)
rbind(sfMmus,sfERCC)
nCountsMmus<-t(t(countsMmus)/sfMmus)
nCountsERCC<-t(t(countsERCC)/sfERCC)
# calculate the sample moments
meansERCC <- rowMeans( nCountsERCC )
varsERCC <- rowVars( nCountsERCC )
cv2ERCC <- varsERCC / meansERCC^2
meansMmus <- rowMeans( nCountsMmus )
varsMmus <- rowVars( nCountsMmus )
cv2Mmus <- varsMmus / meansMmus^2
# fit technical noise
minMeanForFit <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .8 ) )
useForFit <- meansERCC >= minMeanForFit
minMeanForFit
table( useForFit )
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFit] ),
                    cv2ERCC[useForFit] )
fit$coefficients
a0<-unname(fit$coefficients["a0"])
a1<-unname(fit$coefficients["a1tilde"])
# Test for high variance 
minBiolDisp <- .5^2
xi <- mean( 1 / sfERCC )
m <- ncol(countsMmus)
psia1theta <- mean( 1 / sfERCC ) +
  ( coefficients(fit)["a1tilde"] - xi ) * mean( sfERCC / sfMmus )
cv2thA <- coefficients(fit)["a0"] + minBiolDisp + coefficients(fit)["a0"] * minBiolDisp
testDenom <- ( meansMmus * psia1theta + meansMmus^2 * cv2th ) / ( 1 + cv2th/m )
p <- 1 - pchisq( varsMmus * (m-1) / testDenom, m-1 )
# adjust p with the Benjamini-Hochberg method,cut at 10%
padj <- p.adjust( p, "BH" )
sig<-padj<.1
sig[is.na(sig)]<-FALSE
table(sig)
# output table of highly variable genes
log2RelExprMmus<-log2(nCountsMmus/meansMmus)
highVarTable<-data.frame(row.names=NULL,geneID=rownames(countsMmus)[sig],
                         meanNormCount=meansMmus[sig],
                         strongest=factor(colnames(log2RelExprMmus)[
                           apply(log2RelExprMmus[sig,],1,which.max)]),
                         log2RelExprMmus[sig,],
                         check.names=FALSE)
highVarTable[1:10,1:5]
write.csv(highVarTable,file="highly_variant_genes_Mmus.csv",row.names = FALSE)



# Make a plot, highlighting the highly variable ones
colMmus<- "#00207040"
colERCC <- "#70500040"
colMmus <- "#B0901040"
plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 8 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) )
axis( 2, 10^(-2:1), c( "0.01", "0.1", "1", "10" ), las=2 )
abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# Plot the genes, use a different color if they are highly variable
points( meansMmus, cv2Mmus, pch=20, cex=.2,
        col = ifelse( padj < .1, "#C0007090", colMmus ) )
# Add the technical noise fit
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
# Add a curve showing the expectation for the chosen biological CV^2 thershold
lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="#C0007090", lwd=3 )








