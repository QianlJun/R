install.packages("DESeq")
source("https://bioconductor.org/biocLite.R")
biocLite(c("genefilter","EBImage","statmod","topGO","org.At.tair.db"))

library(DESeq)
library(genefilter)
library(EBImage)
library(statmod)
library(topGO)
library(org.At.tair.db)
#options(max.print=300£¬width=100)
sessionInfo()
fullCountTable<-read.csv("C:/Users/Administrator/Desktop/Single-cell/Supplementary_Table_8.csv",header = TRUE,row.names=1)
head(fullCountTable)
countsAll <- fullCountTable[, substr( colnames(fullCountTable), 1, 3 ) == "GL2" ]
colnames(countsAll) <- gsub( "\\.", "-", colnames(countsAll) )
geneTypes <- factor( c( AT="At", pG="pGIBS", EN="HeLa", ER="ERCC" )[
  substr( rownames(countsAll), 1, 2 ) ] )
countsHeLa <- countsAll[ which( geneTypes=="HeLa" ), ]
countsAt <- countsAll[ which( geneTypes=="At" ), ]
countsSp <- countsAll[ which( geneTypes %in% c( "pGIBS", "ERCC" ) ), ]
head(countsHeLa)
head( countsAt )
geneSymbolsAt <- rownames(countsAt)
#rownames(countsAt)
names( geneSymbolsAt ) <- rownames(countsAt)
hasSymbol <- rownames(countsAt) %in% Lkeys( org.At.tairSYMBOL )
symtbl <- toTable( org.At.tairSYMBOL[ rownames(countsAt)[ hasSymbol ] ] )
symtbl <- symtbl[ !duplicated( symtbl$gene_id ), ]
geneSymbolsAt[ symtbl$gene_id ] <- symtbl$symbol
head(geneSymbolsAt)
sfHeLa <- estimateSizeFactorsForMatrix( countsHeLa )
# sfHeLa
sfAt <- estimateSizeFactorsForMatrix( countsAt)
rbind( HeLa = sfHeLa, At = sfAt, ratio = sfAt / sfHeLa )
nCountsHeLa <- t( t(countsHeLa) / sfHeLa )
# t(countsHeLa)
# head(countsHeLa)
# head(nCountsHeLa)
# help(t)
# a<-matrix(1:30,5,6)
# t(a)
nCountsAt <- t( t(countsAt) / sfAt )
colHeLa <- "#00207040"
colAt <- "#70500040"
colAtHi <- "#B0901040"
# pairs( log10( .1 + rbind( nCountsAt, nCountsHeLa ) ), pch=19, cex=.2,
#        col = c( rep( colAt, nrow(nCountsAt) ), rep( colHeLa, nrow(nCountsHeLa) ) ) )
# geneScatterplot( nCountsHeLa[,1], nCountsHeLa[,3],
#                  "normalized read count, GL2 cell 1", "normalized read count, GL2 cell 3",
#                  colHeLa )
meansHeLa <- rowMeans( nCountsHeLa )
#meansHeLa
varsHeLa <- rowVars( nCountsHeLa )
cv2HeLa <- varsHeLa / meansHeLa^2
head(cv2HeLa)
head(varsHeLa)
head(meansHeLa)
#cv2HeLa <- (varsHeLa / meansHeLa)^2
minMeanForFit <- unname( quantile( meansHeLa[ which( cv2HeLa > .3 ) ], .95 ) )
minMeanForFit
# help(unname)
# quantile( meansHeLa[ which( cv2HeLa > .3 ) ], .95 )
# quantile(meansHeLa,0.95)
# summary(meansHeLa)
# length(names(meansHeLa))
useForFit <- meansHeLa >= minMeanForFit
# head(useForFit)
# head(meansHeLa[useForFit])
head(cbind(a0=1,a1tilde=1/meansHeLa[useForFit]))
# a=1
# b=2
# c<-a>=b
# c
# help("glmgam.fit")
# head(cv2HeLa[useForFit])
fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansHeLa[useForFit] ),
                   cv2HeLa[useForFit] )
fit$coefficients
# summary(fit)
xi <- mean( 1 / sfHeLa )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"] - xi )
c( a0, a1 )
# Prepare the plot (scales, grid, labels, etc.)
plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 8 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) )
axis( 2, 10^(-2:1), c( "0.01", "0.1", "1", "10" ), las=2 )
abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# Add the data points
points( meansHeLa, cv2HeLa, pch=20, cex=.2, col=colHeLa )
# Plot the fitted curve
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
# Plot quantile lines around the fit
df <- ncol(countsAt) - 1
lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .975, df ) / df,
       col="#FF000080", lwd=2, lty="dashed" )
lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .025, df ) / df,
       col="#FF000080", lwd=2, lty="dashed" )
meansAt <- rowMeans( nCountsAt )
varsAt <- rowVars( nCountsAt )
cv2At <- varsAt / meansAt^2
psia1theta <- mean( 1 / sfAt ) + a1 * mean( sfHeLa / sfAt )
minBiolDisp <- .5^2
# minBiolDisp
m <- ncol(countsAt)
cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
testDenom <- ( meansAt * psia1theta + meansAt^2 * cv2th ) / ( 1 + cv2th/m )
p <- 1 - pchisq( varsAt * (m-1) / testDenom, m-1 )
padj <- p.adjust( p, "BH" )
padj
sig <- padj < .1
sig[is.na(sig)] <- FALSE
table( sig )
plot( NULL, xaxt="n", yaxt="n",
      log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 8 ),
      xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                       expression(10^4), expression(10^5) ) )
axis( 2, 10^(-2:1), c( "0.01", "0.1", "1", "10" ), las=2 )
abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
# Plot the plant genes, use a different color if they are highly variable
points( meansAt, cv2At, pch=20, cex=.2,
        col = ifelse( padj < .1, "#C0007090", colAt ) )
# Add the technical noise fit, as before
xg <- 10^seq( -2, 6, length.out=1000 )
lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="#C0007090", lwd=3 )



# Analysis of the mouse cells
dataMouse<-read.csv("C:/Users/Administrator/Desktop/Single-cell/Supplementary_Table_5.csv",row.names=1)
dataMouse[ 1:10, 1:5 ]
geneTypes <- factor( c( Mm="Mmus", ER="ERCC" )[
  substr( rownames(dataMouse), 1, 2 ) ] )
countsMmus <- dataMouse[ which( geneTypes=="Mmus" ), -1 ]
countsERCC <- dataMouse[ which( geneTypes=="ERCC" ), -1 ]
lengthsMmus <- dataMouse[ which( geneTypes=="Mmus" ), 1 ]
lengthsERCC <- dataMouse[ which( geneTypes=="ERCC" ), 1 ]
# head(countsMmus)
# head(lengthsMmus)
# head(countsERCC)
# head(lengthsERCC)
sfMmus <- estimateSizeFactorsForMatrix( countsMmus )
sfERCC <- estimateSizeFactorsForMatrix( countsERCC[1:3,] )
sfERCC
rbind( sfMmus, sfERCC )

nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsMmus <- t( t(countsMmus) / sfMmus )
meansERCC <- rowMeans( nCountsERCC )
#meansERCC
varsERCC <- rowVars( nCountsERCC )
cv2ERCC <- varsERCC / meansERCC^2
meansMmus <- rowMeans( nCountsMmus )
varsMmus <- rowVars( nCountsMmus )
cv2Mmus <- varsMmus / meansMmus^2
cv2Mmus
meansERCCPK <- meansERCC / lengthsERCC * 1e3
meansMmusPK <- meansMmus / lengthsMmus * 1e3

minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .8 ) )
useForFitA <- meansERCC >= minMeanForFitA
minMeanForFitA
table( useForFitA )

minMeanForFitB <- unname( quantile( meansERCCPK[ which( cv2ERCC > .3 ) ], .8 ) )
useForFitB <- meansERCCPK >= minMeanForFitB
minMeanForFitB
table( A=useForFitA, B=useForFitB )

fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
                    cv2ERCC[useForFitA] )
fitB <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCCPK[useForFitB] ),
                    cv2ERCC[useForFitB] )
residualA <- var( log( fitted.values(fitA) ) - log( cv2ERCC[useForFitA] ) )
totalA <- var( log( cv2ERCC[useForFitA] ) )
residualB <- var( log( fitted.values(fitB) ) - log( cv2ERCC[useForFitB] ) )
totalB <- var( log( cv2ERCC[useForFitB] ) )
# explained variances of log CV^2 values
c( A = 1 - residualA / totalA,
   B = 1 - residualB / totalB )

minBiolDisp <- .5^2
xi <- mean( 1 / sfERCC )
m <- ncol(countsMmus)
psia1thetaA <- mean( 1 / sfERCC ) +
  ( coefficients(fitA)["a1tilde"] - xi ) * mean( sfERCC / sfMmus )
cv2thA <- coefficients(fitA)["a0"] + minBiolDisp + coefficients(fitA)["a0"] * minBiolDisp
cv2thA
testDenomA <- ( meansMmus * psia1thetaA + meansMmus^2 * cv2thA ) / ( 1 + cv2thA/m )
pA <- 1 - pchisq( varsMmus * (m-1) / testDenomA, m-1 )
padjA <- p.adjust( pA, "BH" )
table( padjA < .1 )
totalVariance <- var( log( cv2ERCC[useForFitA] ) )
m<-ncol(countsERCC)
2/(m-1) / totalVariance

