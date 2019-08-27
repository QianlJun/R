##p299
library(psych)
fa.parallel(USJudgeRatings[,-1],fa="pc",n.iter=100,
            show.legend = FALSE,main="Scree plot with parallel analysis")

##p300
library(psych)
pc<-principal(USJudgeRatings[,-1],nfactors = 1)
pc

##p301
library(psych)
fa.parallel(Harman23.cor$cov,n.obs = 302,fa="pc",n.iter=100,
            show.legend = FALSE,main="Scree plot with parallel analysis")

##p302
library(psych)
pc<-principal(Harman23.cor$cov,nfactors = 2,rotate = "none")
pc

##p303
rc<-principal(Harman23.cor$cov,nfactors = 2,rotate = "varimax")
rc

##p304
library(psych)
pc<-principal(USJudgeRatings[,-1],nfactors = 1,score=TRUE)
head(pc$scores)

##p306
options(digits = 2)
covariances<-ability.cov$cov
correlations<-cov2cor(covariances)
correlations

library(psych)
covariances<-ability.cov$cov
correlations<-cov2cor(covariances)
fa.parallel(correlations,n.obs=112,fa="both",n.iter=100,
            main="Scree plots with parallel analysis")

##p308
fa<-fa(correlations,nfactors=2,rotate="none",fm="pa")
fa

fa.varimax<-fa(correlations,nfactors=2,rotate="varimax",fm="pa")
fa.varimax

##P309
fa.promax<-fa(correlations,nfactors=2,rotate="promax",fm="pa")
fa.promax

##p310
fsm<-function(oblique){
  if (class(oblique)[2]=="fa" & is.null(oblique$Phi)){
    warning("Object doesn't look like oblique EFA")
  }
  else {
    P<-unclass(oblique$loading)
    F<-P %*% oblique$Phi
    colnames(F)<-c("PA1","PA2")
    return(F)
  }
}
fsm(fa.promax)

factor.plot(fa.promax,labels=rownames(fa.promax$loadings))
fa.diagram(fa.promax,simpe=FALSE)
