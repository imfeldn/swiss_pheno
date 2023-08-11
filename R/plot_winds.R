

### wind arrow plot
quiver<- function(u,v,x.val,y.val,scale=scale,length=leng,col,lwd=1)
{
  xx <- matrix(ncol=length(x.val),rep((x.val),length(y.val)),nrow=length  (y.val),byrow=TRUE)
  yy <- matrix(ncol=length(x.val),rep((y.val),length(x.val)),nrow=length    (y.val),byrow=FALSE)
  # I removed scaling from this function...
  speed <- sqrt(u*u+v*v)
  maxspeed <- max(speed)
  u <- u *scale/maxspeed
  v <- v *scale/maxspeed
  arrows(xx,yy[,(length(yy[1,])):1],xx+t(u),yy[,(length(yy[1,])):1]+t(v),length=length*min(par.uin()),lwd=lwd, col =col)
}

par.uin <- function()# determine scale of inches/userunits in x and y
{
  u <- par("usr")
  p <- par("pin")
  c(p[1]/(u[2] - u[1]), p[2]/(u[4] - u[3]))
}

quiver.leg <- function(u,v,x.val,y.val,scale=scale,length=leng){
  xx <- matrix(ncol=length(x.val),rep((x.val),length(y.val)),nrow=length	(y.val),byrow=TRUE)
  yy <- matrix(ncol=length(x.val),rep((y.val),length(x.val)),nrow=length		(y.val),byrow=FALSE)
  speed <- sqrt(u*u+v*v)
  maxspeed <- max(speed)
  u <- u*scale/maxspeed
  v <- v*scale/maxspeed
  arrows(xx,yy[,length(yy[,1]):1],xx+u,yy[,length(yy[,1]):1]+v,length=leng*min(par.uin()),lwd=0.8)
}
