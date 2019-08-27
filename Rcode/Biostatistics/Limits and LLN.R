nosim<-1000
x<-cumsum(rnorm(nosim))/1:nosim
plot(1:nosim,x,type="l",xlab="Iteration",ylab="Average",frame=FALSE)
abline(h=0,lty=2)