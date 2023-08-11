

# function to make colors transparent
t_col<-function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# function to create a customized image scale
image_scale <- function(breaks,col, axis.pos=1, add.axis=TRUE, xlab=NULL,ylab=NULL,scale_lab, key.extend=c(FALSE,FALSE),useraxis=NULL, equidist=FALSE,mar=c(5,1,5,3),las=1,...){

  poly <- vector(mode="list", length(col))
  for(i in seq(poly)[1:length(poly)]){
    if(equidist==TRUE){
      poly[[i]] <- c(i+key.extend[1],i+key.extend[1]+1, i+key.extend[1]+1, i+key.extend[1])
      #poly[[i]] <- c(i,i+1, i+1, i)
    } else  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }

  ifelse(equidist,height<-c((key.extend[1]+1):(length(breaks)-key.extend[2])),height<-breaks)
  #ifelse(equidist,height<-c((key.extend[1]+1):(length(breaks)-key.extend[2])),height<-breaks[(key.extend[1]+1):(length(breaks)-key.extend[2])])
  if(axis.pos %in% c(1,3)){ylim<-c(0,1);xlim=range(height) }
  if(axis.pos %in% c(2,4)){ylim<-range(height); xlim<-c(0,1)}
  if(missing(ylab)){ylab=""}
  if(missing(xlab)){xlab=""}
  par(mar=mar)
  plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab=xlab, ylab=ylab, xaxs="i", yaxs="i") #, ...)
  for(i in seq(poly)){
    if(axis.pos %in% c(1,3)){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(axis.pos %in% c(2,4)){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }

  if(axis.pos %in% c(2,4) & is.null(useraxis)){
    mtext(side=4,at=breaks[round(length(breaks)/2)],scale_lab,line=1.5,cex=1,las=las, padj = 1)
  }
  if (any(key.extend)) {
    dl <- diff(breaks[(length(breaks)-1):length(breaks)])
    last <- length(breaks)
    xi <- 0
    xa <- 1
    apex <- 1
    clipmax <- apex + (0.05*apex)
    if(axis.pos %in% c(1,3)){
      # rect(breaks[(1+key.extend[1]):(last-1-key.extend[2])],xi,
      #   breaks[(2+key.extend[1]):(last-key.extend[2])],xa,
      #   col = col[(1+key.extend[1]):(length(col)-key.extend[2])])
      clip(breaks[1]-(dl*clipmax), breaks[last]+(dl*clipmax),xi,xa )
    } else {
      # rect(xi, breaks[(1+key.extend[1]):(last-1-key.extend[2])],
      # xa, breaks[(2+key.extend[1]):(last-key.extend[2])],
      # col = col[(1+key.extend[1]):(length(col)-key.extend[2])])
      clip(xi,xa, breaks[1]-(dl*clipmax), breaks[last]+(dl*clipmax))
    }

    if (key.extend[1]) {
      lower_break <- breaks[1]-dl
      if(axis.pos %in% c(1,3)){
        polygon( c(breaks[1], breaks[1],
                   lower_break-(dl*apex)),c(xi,xa,xa/2),
                 col = col[1], border="black")
      } else {polygon(c(xi,xa,xa/2),
                      c(breaks[1], breaks[1], lower_break),
                      col = col[1], border=NA)}
    }
    if (key.extend[2]) {
      upper_break <- breaks[last]+dl
      if(axis.pos %in% c(1,3)){
        polygon( c(breaks[last]+(dl), breaks[last], breaks[last],
                   breaks[last]+(dl), upper_break+(dl*apex)),
                 c(xi,xi,xa,xa,xa/2),
                 col = col[length(col)], border=NA)

      } else {polygon(c(xi,xa/2,xa),
                      c(breaks[last], upper_break, breaks[last]),
                      col = col[last-1], border=NA)}
    }
    new_low<- ifelse(key.extend[1], lower_break, breaks[1])
    new_up<- ifelse(key.extend[2], upper_break, breaks[last])

    if(axis.pos %in% c(1,3)){
      polygon(c(breaks[last],new_up, breaks[last],breaks[1], new_low, breaks[1]),c(xi,xa/2,xa,xa,xa/2,xi),
              col = NA, border="black" )
    } else{
      polygon(c(xi,xa/2,xa,xa,xa/2,xi),c(breaks[last],new_up, breaks[last],breaks[1], new_low, breaks[1]),
              col = NA, border="black" )


    }

  } else
  box()
  if(add.axis) axis(axis.pos)
  if(!is.null(useraxis)){

    axis(side=useraxis$axis.pos,at=useraxis$ats, labels=useraxis$labs, line=0, cex.axis=useraxis$cex,las=1)
    #axis(side=useraxis$axis.pos,at=useraxis$ats, labels=useraxis$labs, line=0, cex.lab=useraxis$cex,las=1)
    mtext(side=useraxis$axis.pos2,at=breaks[round(length(breaks)/2)],scale_lab,line=useraxis$line,cex=useraxis$cex, adj = useraxis$adj)

  }

}
