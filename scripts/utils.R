short.abline<-function(x,y,col="black",lwd=1)
{
   lines(x,y,col=col,lwd=lwd)
}


double.lines<-function(x,y,bcol="black",fcol,lwd,lty="solid",pretty=TRUE)
{

  if (pretty)
  {
    lines(x,y,col=bcol,lwd=lwd+2)
    lines(x,y,col=fcol,lwd=lwd,lty=lty)
  }
  else
    lines(x,y,col=fcol,lwd=2)
}

nice.ts.plot<-function(x,y,label=NULL,legend=NULL,bcol="blue",fcol="lightblue",lwd=5,ylab=NULL)
{
   nlines <- 1
   if (!is.null(ncol(y)))
     nlines = ncol(y)

   yrange <- c(0,1.2*max(y,na.rm=TRUE))
#  yrange <- c(min(y,na.rm=TRUE),1.2*max(y,na.rm=TRUE))
   q <- quantile(y,probs=c(0.01,0.99),na.rm=TRUE)
#  yrange <- c(q[1],(1.6*q[2]))
#  print(yrange)
   ytic <- pretty(c(yrange[1],yrange[2]),5)
#  ytic <- seq(yrange[1],yrange[2],((yrange[2]-yrange[1])/5))
#  print(ytic)
   xrange <- c(min(x,na.rm=TRUE),max(x,na.rm=TRUE))
#  print(xrange)
   if (!is.null(legend))
   {
   # print(paste(legend,nchar(legend)))
     xrange[2] = xrange[2]+nchar(legend)*1.5
   # print(xrange)
   }

   old.par <- par(no.readonly = TRUE) 
   par(mar=c(3,4,1,1)+0.1)
   plot(xrange,yrange,type='n',axes=FALSE,ann=FALSE,yaxs='i')
   for (i in 2:length(ytic))
   {
      short.abline(c(par("usr")[1],range(x)[2]),c(ytic[i],ytic[i]),col="lightgray",lwd=2)
   #  abline(h=ytic[i],col="lightgray",lwd=2)
   } 
   if (nlines > 1)
   {
      for (n in 1:nlines)
      {
        double.lines(x,y[,n],bcol=bcol,fcol=fcol,lwd=lwd)
      }
   } 
   else
   {
      double.lines(x,y,bcol=bcol,fcol=fcol,lwd=lwd)
   }
   if (!is.null(legend))
      text(x[length(x)],y[length(y)],legend,col=bcol,pos=4,offset=0.25,cex=1.0)
   if (!is.null(label))
      title(main=label,line=-0.75) #outer=FALSE)

   axis(1,lwd=0)
   abline(h=par("usr")[3],lwd=3)

   axis(2,lwd=0,las=1,at=ytic[2:length(ytic)])
   abline(v=par("usr")[1],lwd=3)
   mtext(ylab,side=2,line=3)
#  par(old.par)
}

save.png.plot<-function(root,width=6.5,height=4.5)
{
  graphics.root <-paste("../Reports/graphics/",root,sep="")
  file.png <-paste(graphics.root,".png",sep="")
  file.pdf <-paste(graphics.root,".pdf",sep="")
  dev.copy2pdf(file=file.pdf,width=width,height=height)
# pdfcrop --margins 0 file.pdf
# file.eps <-paste(graphics.root,".eps",sep="")
# dev.copy2eps(file=file.eps,width=width,height=height)
  cmd <- paste("convert -antialias -density 300",file.pdf,file.png,sep=" ")
  system(cmd)
  print(paste("Plot saved as ",file.pdf," and converted to ", file.png,sep=""),
              quote=FALSE)
# system(paste("rm -fv",file.pdf))
}

