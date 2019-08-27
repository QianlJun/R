library(rgl)  
#预测脚本  
predictgrid<-function(model,xvar,yvar,zvar,res=16,type=NULL){  
  xrange<-range(model$model[[xvar]])  
  yrange<-range(model$model[[yvar]])  
  newdata<-expand.grid(x=seq(xrange[1],xrange[2],length.out=res),  
                       y=seq(yrange[1],yrange[2],length.out=res))  
  names(newdata)<-c(xvar,yvar)  
  newdata[[zvar]]<-predict(model,newdata=newdata,type=type)  
  newdata  
}  
#x,y,z转为列表  
df2mat<-function(p,xvar=NULL,yvar=NULL,zvar=NULL){  
  if(is.null(xvar)) xvar<-names(p)[1]  
  if(is.null(yvar)) yvar<-names(p)[2]  
  if(is.null(zvar)) zvar<-names(p)[3]  
  x<-unique(p[[xvar]])  
  y<-unique(p[[yvar]])  
  z<-matrix(p[[zvar]],nrow=length(y),ncol=length(x))  
  m<-list(x,y,z)  
  names(m)<-c(xvar,yvar,zvar)  
  m  
}  
#交错出现两个向量元素  
interleave<-function(v1,v2) as.vector(rbind(v1,v2))  

m<-mtcars  
mod<-lm(mpg~wt+disp+wt:disp,data=m)  
m$pred_mpg<-predict(mod)  
mpgrid_df<-predictgrid(mod,'wt','disp','mpg')  
mpgrid_list<-df2mat(mpgrid_df)  

plot3d(mtcars$wt,mtcars$disp,mtcars$mpg,xlab='',ylab='',zlab='',axes=FALSE,size=.5,type='s',lit=FALSE)  

spheres3d(m$wt,m$disp,m$pred_mpg,alpha=0.4,type='s',size=0.5,lit=FALSE)  

segments3d(interleave(m$wt, m$wt),  
           interleave(m$disp, m$disp),  
           interleave(m$mpg, m$pred_mpg),  
           alpha=0.4,col='red'  
)  
#预测曲面  
surface3d(mpgrid_list$wt,mpgrid_list$disp,mpgrid_list$mpg,alpha=.4,front='lines',back='lines')  
#其他设置  
rgl.bbox(color='grey50',emission='grey50',xlen=0,ylen=0,zlen=0)  
rgl.material(color='black')  
axes3d(edges=c('x--','y+-','z--'), cex=.75)  
mtext3d('Society',edge='x--',line=2)  
mtext3d('Pollution',edge='y+-',line=3)  
mtext3d('Economy',edge='z--',line=3)  