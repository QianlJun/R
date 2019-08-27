##p265_1
library(coin)
score<-c(40,57,45,55,58,57,64,55,62,65)
treatment<-factor(c(rep("A",5),rep("B",5)))
mydata<-data.frame(treatment,score)
t.test(score~treatment,data=mydata,var.equal=TRUE)
oneway_test(score~treatment,data=mydata,distribution="exact")

##p265_3
library(multcomp)
set.seed(1234)
oneway_test(response~trt,data=cholesterol,
            distribution=approximate(B=9999))

##p274
rsq<-function(formula,data,indices){
  d<-data[indices,]
  fit<-lm(formula,data=d)
  return(summary(fit)$r.square)
}
library(boot)
set.seed(1234)
results<-boot(data=mtcars,statistic=rsq,
              R=1000,formula=mpg~wt+disp)
print(results)
plot(results)
boot.ci(results,type=c("perc","bca"))
