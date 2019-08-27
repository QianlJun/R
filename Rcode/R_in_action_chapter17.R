pkgs<-c("rpart","rpart.plot","party","randomForest","e1071")
install.packages(pkgs,depend=TRUE)

##p361
loc<-"http://archive.ics.uci.edu/ml/machine-learning-databases/"
ds<-"breast-cancer-wisconsin/breast-cancer-wisconsin.data"
url<-paste(loc,ds,sep="")
breast<-read.table(url,sep=",",header=FALSE,na.strings="?")
names(breast)<-c("ID","clumpThickness","sizeUniformity",
                 "shapeUniformity","maginalAdhesion",
                 "singleEpitheliaCellSize","bareNuclei",
                 "blandChromatin","normalNucleoli","mitosis","class")
df<-breast[-1]
df$class<-factor(df$class,levels=c(2,4),
                 labels=c("benign","malignant"))
set.seed(1234)
train<-sample(nrow(df),0.7*nrow(df))
df.train<-df[train,]
df.validate<-df[-train,]
table(df.train$class)
table(df.validate$class)

##p362
fit.logit<-glm(class~.,data=df.train,family=binomial())
summary(fit.logit)
prob<-predict(fit.logit,df.validate,type="response")
logit.pred<-factor(prob>.5,levels=c(FALSE,TRUE),
                     labels=c("benign","malignant"))
logit.perf<-table(df.validate$class,logit.pred,
                  dnn=c("Actual","Predicted"))
logit.perf

##p364
library(rpart)
set.seed(1234)
dtree<-rpart(class~.,data=df.train,method="class",
             parms=list(split="information"))
dtree$cptable
plotcp(dtree)
dtree.pruned<-prune(dtree,cp=.0125)
library(rpart.plot)
prp(dtree.pruned,type=2,extra=104,
    fallen.leaves=TRUE,main="Decision Tree")
dtree.pred<-predict(dtree.pruned,df.validate,type="class")
dtree.perf<-table(df.validate$class,dtree.pred,
                  dnn=c("Actual","Predicted"))
dtree.perf

##p367
library(party)
fit.ctree<-ctree(class~.,data=df.train)
plot(fit.ctree,main="Conditional Inference Tree")
ctree.pred<-predict(fit.ctree,df.validate,type="response")
ctree.perf<-table(df.validate$class,ctree.pred,
                  dnn=c("Actual","Predicted"))
ctree.perf

##p369
library(randomForest)
set.seed(1234)
fit.forest<-randomForest(class~.,data=df.train,
                         na.action=na.roughfix,
                         importance=TRUE)
fit.forest
importance(fit.forest,type=2)
forest.pred<-predict(fit.forest,df.validate)
forest.perf<-table(df.validate$class,forest.pred,
                   dnn=c("Actual","Predicted"))
forest.perf

##p372
library(e1071)
set.seed(1234)
fit.svm<-svm(class~.,data=df.train)
fit.svm
svm.pred<-predict(fit.svm,na.omit(df.validate))
svm.perf<-table(na.omit(df.validate)$class,
                        svm.pred,dnn=c("Actual","Predicted"))
svm.perf

##p373
set.seed(1234)
tuned<-tune.svm(class~.,data=df.train,
                gamma=10^(-6:1),
                cost-10^(-10:10))
tuned
fit.svm<-svm(class~.,data=df.train,gamma=0.01,cost=1)
svm.pred<-predict(fit.svm,na.omit(df.validate))
svm.perf<-table(na.omit(df.validate)$class,
                svm.pred,dnn=c("Actual","Predicted"))
svm.perf

##p374
performance<-function(table,n=2){
  if(!all(dim(table)==c(2,2)))
    stop("Must be a 2x2 table")
  tn=table[1,1]
  fp=table[1,2]
  fn=table[2,1]
  tp=table[2,2]
  sensitivity=tp/(tp+fn)
  specificity=tn/(tn+fp)
  ppp=tp/(tp+fp)
  npp=tn/(tn+fn)
  hitrate=(tp+tn)/(tp+tn+fp+fn)
  result<-paste("Sensitivity= ",round(sensitivity,n),
                "\nSpecificity = ",round(specificity,n),
                "\nPositive Predictive Value = ",round(ppp,n),
                "\nNegative Predictive Value = ",round(npp,n),
                "\nAccuracy = ",round(hitrate,n),"\n",sep=" ")
  cat(result)
}

##p375
performance(logit.perf)
performance(dtree.perf)
performance(ctree.perf)
performance(forest.perf)
performance(svm.perf)

##p377
library(rattle)
rattle()
loc<-"http://archive.ics.uci.edu/ml/machine-learning-databases/"
ds<-"pima-indians-diabetes/pima-indians-diabetes.data"
url<-paste(loc,ds,sep="")
diabetes<-read.table("F:/workdata/pima-indians-diabetes.data.csv", sep=",",header=FALSE)
names(diabetes)<-c("npregant","plasma","bp","triceps",
                   "insulin","bmi","pedigree","age","class")
diabetes$class<-factor(diabetes$class,levels=c(0,1),
                       labels=c("normal","diabetic"))
library(rattle)
rattle()

