library(rgl)
library(mvtnorm)
sigma<-matrix(c(8,0.25*8,0.25*8,8),2,2)
sigma
xvals<-seq(-10,10,length=100)
yvals<-seq(-10,10,length=100)
zvals<-apply(expand.grid(xvals,yvals),1,
             function(w) dmvnorm(w,mean = c(0,0),sigma=sigma))
persp3d(x=xvals,y=yvals,z=zvals,col="lightblue")
planes3d(0,1,0,-5,col=grey(0.8))
curve(x^3-x^4,0,1)
curve(x^2,-2,2)
curve(dnorm(x,mean=0.712,sd=(0.031)^0.5),from = -1,to=2)
qnorm(0.95)
