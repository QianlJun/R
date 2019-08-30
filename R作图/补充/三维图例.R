library(lattice)
x=seq(-pi,pi,len=20)
y=x
g=expand.grid(x=x,y=y)
g$z=sin(sqrt(g$x^2+g$y^2))
wireframe(z~x*y,data=g,drape=T,aspect=c(3,1),colorkey=T,main=expression(z=sin(sqrt(x^2+y^2))))