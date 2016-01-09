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

nice.ts.plot<-function(x,y,label=NULL,legend=NULL,bcol="blue",fcol="lightblue",
              lwd=5, ylab=NULL,xlab=NULL,ylim=NULL)
{
   nlines <- 1
   if (!is.null(ncol(y)))
     nlines = ncol(y)

   if (is.null(ylim))
      yrange <- c(0,1.2*max(y,na.rm=TRUE))
   else
      yrange=ylim
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
   #  print(xrange)
   #  print(paste(legend,max(nchar(legend))))
      xrange[2] = xrange[2]+max(nchar(legend))*1.2
   #  print(xrange)
   }

   old.par <- par(no.readonly = TRUE) 
 # par(mar=c(3,4,2,1)+0.1)
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
   {
      if (nlines > 1)
      {
         for (n in 1:nlines)
         {
            text(x[length(x)],y[nrow(y),n],legend[n],col=bcol,pos=4,offset=0.25,cex=1.0)
         }
      }
      else
        text(x[length(x)],y[length(y)],legend,col=bcol,pos=4,offset=0.25,cex=1.0)
   }
   if (!is.null(label))
      title(main=label,line=-0.75) #outer=FALSE)

   axis(1,lwd=0)
   abline(h=par("usr")[3],lwd=3)
   mtext(xlab,side=1,line=2.5)

   axis(2,lwd=0,las=1,at=ytic[2:length(ytic)])
   abline(v=par("usr")[1],lwd=3)
   mtext(ylab,side=2,line=3)

   return(xrange)
#  par(old.par)
}


label.panel=function(label)
{
   legend("topleft",bty='n',cex=1.6,legend=label,text.font=2)
}

sd.bars = function(x,y,s,col="orange",lwd=3,lty="dashed")
{
#  print("sd.bars:")
#  print(x)
#  print(y)
#  print(s)
   np = length(x)
   for (p in 1:np)
   {
      lines(c(x[p],x[p]),c((y[p]-s[p]),(y[p]+s[p])),
       col=col,lwd=lwd,lty=lty)
   #  double.lines(c(x[p],x[p]),c((y[p]-s[p]),(y[p]+s[p])),
   #   bcol="blue",fcol="lightblue",lwd=lwd)
   }
}


save.png.plot<-function(root,width=6.5,height=4.5)
{
# graphics.root <-paste("../Reports/graphics/",root,sep="")
  graphics.root <-paste("./",root,sep="")
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

plot.error=function(x,y,sd,bcol,fcol,mult=2)
{
   if (capabilities("cairo"))
   {
      sdyu = exp(log(y)+mult*sd)
      sdyl = exp(log(y)-mult*sd)
      frgb = col2rgb(fcol)/255
      polygon(c(x,rev(x)),c(sdyl,rev(sdyu)),
              border=bcol,lty="dashed",lwd=1,
              col=rgb(frgb[1],frgb[2],frgb[3],0.5))
   }
   else
      polygon(c(x,rev(x)),c(sdyl,rev(sdyu)),
              border=bcol,lty="dashed",lwd=1,col=fcol)

}


show.block.number=function(block.number,x,line=3,col="black")
{
   mtext(text=paste("(",block.number,")",sep=""),side=1,line=line,
          at=c(x,0),cex=0.8,col=col)
}

plot.catches=function(t,obs,pred,sd=NULL,block=NULL)
{
#  print(head(t))
#  print(head(obs))
#  print(head(pred))
   ntime = length(t)
   ngear = ncol(obs)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(3,4.5,0,4)+0.1)
   np = ngear+1
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   lwd = 3
   title.line = -1
   sum.obs = vector(length=ntime,mode="numeric")
   sum.pred = vector(length=ntime,mode="numeric")
   for (g in 1:ngear)
   {
      nice.ts.plot(t,pred[,g],
         bcol="darkgreen",fcol="lightgreen",lwd=lwd,ylab="Catch (mt)")
      if (!is.null(sd))
        plot.error(t,pred[,g],sd,
                 bcol="darkgreen",fcol="lightgreen")
      lines(t,pred[,g],col="darkgreen",lwd=lwd+2)
      points(t,obs[,g],col= "darkgreen",pch=3,cex=3)
      title(main=gear.names[g],line=title.line)

      sum.obs  = sum.obs  + obs[,g]
      sum.pred = sum.pred + pred[,g]
   }

   nice.ts.plot(t,sum.pred,
                 bcol="darkgreen",fcol="lightgreen",lwd=lwd,ylab="Catch (mt)")
   lines(t,sum.pred,col="darkgreen",lwd=lwd+2)
   points(t,sum.obs,col= "darkgreen",pch=3,cex=3)
   title(main="All Fleets",line=title.line)

   if (!is.null(block))
      show.block.number(block,t[1],line=2)

   par(old.par)
}

plot.biomass=function(x,y,K=NULL,sd=NULL,block=NULL,propL=NULL,
                      forcing=NULL,yrange=NULL)
{
#  print(length(x))
#  print(head(x))
#  print(dim(y))
#  print(head(y))
   ntime = length(x)
#  print(ntime)

   legend = c(parse(text=paste("N","[1]")),
              parse(text=paste("N","[2]")),
              parse(text=paste("N","[1]","+N","[2]")))

   options(scipen=6)
   if (is.null(yrange))
      yrange = c(0,1.2*max(y,na.rm=TRUE))
   xrange=nice.ts.plot(x,y,ylab="Biomass (mt)",ylim=yrange)
   print(ncol(y))

   if (!is.null(K))
   {
      lines(x,K,lwd=2,lty="dotdash",col="blue",xlim=xrange)
      text(x[1],K[1],"K ",adj=c(1.5,0.5),col="blue",cex=1.2)
   }

#  sdlogNN = sqrt(4.0*sd*sd) # this is probably not correct
#  print("N3 error")
#  if ( (ncol(y) > 2) && (!is.null(sd)) )
#     plot.error(x,y[,3],sdlogNN, bcol="blue",fcol="lightblue")

   if (!is.null(sd))
   {
      print("N1 error")
      plot.error(x,y[,1],sd, bcol="blue",fcol="lightblue")
   }
   if ( (ncol(y) > 1) && (!is.null(sd)) )
   {
      print("N2 error")
      plot.error(x,y[,2],sd, bcol="blue",fcol="lightblue")
      lines(x,y[,1],col="blue",lwd=5)
   }
   print("N1")
   lines(x,y[,1],col="blue",lwd=5)
   print("N1 label")
   text(x[ntime],y[ntime,1],parse(text=paste("N","[1]")),adj=c(-0.5,0.5),
        col="blue",cex=1.2)

   if (ncol(y) > 1)
   {
      print("N2")
      lines(x,y[,2],col="blue",lwd=5)
#  print("N2 label")
      text(x[ntime],y[ntime,2],parse(text=paste("N","[2]")),adj=c(-0.5,0.5),col="blue",cex=1.2)
   }

   if (ncol(y) > 2)
   {
      print("N3")
      lines(x,y[,3],col="blue",lwd=5)
      text(x[ntime],y[ntime,3],parse(text=paste("N","[1]","+","N","[2]")),adj=c(0.5,0.5),col="blue",cex=1.2)
   }
   if (!is.null(propL))
   { 
      print("propL")
      par("new"=TRUE)
      plot(x,propL,lwd=3,type='l',col="red",ylim=c(0,1),
            ann=FALSE,axes=FALSE,xlim=xrange)
      text(x[ntime],propL[ntime]," p",adj=c(0,0),col="red")
      abline(h=0.9,lwd=2,lty="dotdash",col="red")
      axis(4,col="red",ylab="p",col.axis="red",line=-2.0)
      mtext("p",side=4,col="red",line=-1.0)
   }
   if ( (!is.null(forcing)) && (max(forcing)> 1) )
   {
      print("forcing")
      par("new"=TRUE)
      frange = c(0,1.2*max(forcing,na.rm=TRUE))
      ftic <- pretty(c(frange[1],frange[2]),5)
      plot(xrange,frange,type='n',axes=FALSE,ann=FALSE,yaxs='i')
      lines(x,forcing,lwd=3,col="purple")
   #  plot(x,forcing,lwd=3,type='l',col="purple",
   #       ylim=frange,ann=FALSE,axes=FALSE,xlim=xrange)
      text(x[1],forcing[1],"I ",adj=c(1.5,0.5),col="purple",cex=1.2)
      axis(4,lwd=1,at=ftic[2:length(ftic)],col="purple",col.axis="purple")
      abline(v=par("usr")[2],lwd=2,col="purple")
      mtext("Forcing (I)",side=4,col="purple",line=3)
    }
    if (!is.null(block))
       show.block.number(block,x[1])
}


plot.production=function(Fmort, obsC, predC, t, r, K, block=NULL)
{
    F.max=max(Fmort,na.rm=TRUE)
    Fyield = seq(0,F.max,0.01*F.max)
    yield = Fyield*K*(1.0-Fyield/r) # equilibirum yield at F
    #  obsC = rowSums(dat[,obsC.ndx])
    #  predC = rowSums(dat[,predC.ndx])
    xrange = c(0,F.max)
    yrange = c(0,max(obsC,predC,yield,na.rm=TRUE))
    Flab = parse(text=paste("Total~Fishing~Mortality~(","y^-1",")",sep=""))
    plot(xrange,yrange,type='n', xlab=Flab, ylab="Total Yield (mt)")

    double.lines(Fmort,predC,bcol="darkgreen",fcol="lightgreen",lwd=5) 
    points(Fmort,obsC,col= "darkgreen",pch=3,cex=2)
    points(Fmort,predC,col="darkgreen",pch=16)
   #   wmaxC = which(obsC==max(obsC,na.rm=TRUE))
   #   lines(c(0,Fmort[wmaxC]),c(0,obsC[wmaxC]),col="red",lty="longdash")
    lines(Fyield,yield,col="red",lwd=3,lty="longdash")
    wy5 = which((floor(t%%5))==0)
    text(x=Fmort[wy5],y=obsC[wy5],labels=floor(t[wy5]),
         pos=4,offset=0.5,cex=0.8)
    if (!is.null(block))
       show.block.number(block,Fyield[1],line=4)
}

#read.diagnostics=function(file="xssams_program.log",ntime=61,dt=1,ngear=5,block=NULL)
get.diagnostics=function(log,ntime=61,dt=1,ngear=5,block=NULL,mtype)
{
   get.numeric.field = function(what, text)
   {
      old.opt=getOption("warn")
      options("warn"=-1)
      ndx = grep(what,text)
      lwhat = nchar(what)
      value = 0.0
      if (substr(what,lwhat,lwhat) == ":")
         ndx = ndx + 1
      else
         ndx = ndx + 2

      value = as.numeric(text[ndx])
      wv = which(!is.na(value))
      options("warn"=old.opt)
      return(value[wv])   
   } 



#  print(paste("Scanning file",file))
#  log = scan(file,what="character")
## print(length(log))
   res = grep("Residuals:",log)
   nblock = length(res)
   print(paste(nblock,"residual blocks found"))

   logsdlogF = get.numeric.field("logsdlogF:",log)   
   sdlogF = exp(logsdlogF)

   if (is.null(block))
      block = nblock

   ests = list(logT12=NA, T12=NA, logT21=NA, T21=NA, logr=NA,
               r=NA, r_prior=NA, sdr_prior=NA, logB1=NA, logdB1K=NA,
               B1=NA, K=NA, MSY=NA, Fmsy=NA, logsdlogProc=NA,
               sdlogPop=NA, logsdlogYield=NA, sdlogYield=NA, pcon=NA,
               qProp=NA, logsdlogF=NA, sdlogF=NA, logsdlogPop=NA,
               resid=NA)

#  print(ests)

   ests$logT12=get.numeric.field("^logT12",log)[block]
   ests$T12=get.numeric.field("^T12",log)[block]
   ests$logT21=get.numeric.field("^logT21",log)[block]
   ests$T21=get.numeric.field("^T21",log)[block]
   ests$logr=get.numeric.field("^logr",log)[block]
   ests$r=get.numeric.field("^r",log)[block]
   ests$r_prior=get.numeric.field("^r_prior",log)[block]
   ests$sdr_prior=get.numeric.field("^sdr_prior",log)[block]
   ests$logB1=get.numeric.field("^logB1",log)[block]
   ests$logdB1K=get.numeric.field("^logdB1K",log)[block]
   ests$B1=get.numeric.field("^B1",log)[block]
   ests$K=get.numeric.field("^K",log)[block]
   ests$MSY=get.numeric.field("^MSY",log)[block]
   ests$Fmsy=get.numeric.field("^Fmsy",log)[block]
   ests$logsdlogProc=get.numeric.field("^logsdlogProc:",log)[block]
   ests$sdlogPop=get.numeric.field("^sdlogPop:",log)[block]
   ests$logsdlogYield=get.numeric.field("^logsdlogYield:",log)[block]
   ests$sdlogYield=get.numeric.field("^sdlogYield:",log)[block]
   ests$pcon=get.numeric.field("^pcon",log)[block]
   ests$qProp=get.numeric.field("^qProp",log)[block]
   ests$logsdlogF=get.numeric.field("^logsdlogF:",log)[block]
   ests$sdlogF=get.numeric.field("^sdlogF:",log)[block]
   ests$logsdlogPop=get.numeric.field("^logsdlogPop:",log)[block]


   # npop?
   if (mtype == "x")
      ncol = (3*ngear+6)
   else if (mtype == "i")
      ncol = (3*ngear+4)
   else
      stop(paste("Unknown model type(",mtype,
                 ") passed to get.diagnostics(...)",sep=""))
 
   resid = matrix(nrow=ntime,ncol=ncol)
   cnames = vector(length=ncol)

   fc = res[block]
   for (c in 1:ncol)
   {
      fc = fc + 1
      cnames[c] = log[fc]
   }
   colnames(resid) = cnames
#  print(colnames(resid))
 
   for (t in 1:ntime)
   {
      for (c in 1:ncol)
      {
         fc = fc + 1
         resid[t,c] = as.numeric(log[fc],ngear)
      }
   }

#  print(paste("Block",block,":"))
#  print(ests)
#  print(head(resid))

   ests$resid=as.data.frame(resid)

   return(ests)
 
}


plot.biomass.array=function(path.list=c("../run-issams-dev/issams-dev.rep",
                                         "../run-issams/issams.rep",
                                         "../run-xssams/1/xssams.rep"),
                             ntime=61,dt=1,ngear=5,mtype=c("i","i","x"))
{
   npath = length(path.list)
   print(npath)
   width = 9.0
   height = 11.0
   yrange=c(0,30000)
   start.year = 1952

   old.par = par(no.readonly = TRUE,"new"=FALSE) 

   x11(width=width,height=height,title="Biomass")
   par(mar=c(3,4.5,0,5)+0.1,"new"=FALSE)
   lm = layout(matrix(c(1:npath),ncol=1,byrow=TRUE))
   layout.show(lm)
   biomass.dev = dev.cur()

   x11(width=width,height=height,title="Production",xpos=100)
   par(mar=c(3,4.5,0,0)+0.1,"new"=FALSE,las=1)
   lm = layout(matrix(c(1:npath),ncol=1,byrow=TRUE))
   layout.show(lm)
   prod.dev = dev.cur()

   for (p in 1:npath)
   {
      path=path.list[p]
      print(path)
      pv = unlist(strsplit(path,"[\\./]"))
      model = pv[length(pv)-1]
      print(model)
      if (model == "issams-dev")
         #     parse(text=paste("N","[1]","+","N","[2]"))
         legend = parse(text=paste("B","[1]","~d","~~no~index"))
      else if (model == "issams")
         legend = parse(text=paste("MSY~F","[msy]","~~indexed"))
      else if (model == "xssams")
         legend = parse(text=paste("MSY~F","[msy]","~~immigration"))
      else
         legend = model

      log = scan(path,what="character")
      res = grep("Residuals:",log)
      block = length(res)
      tmp=get.diagnostics(log,ntime=ntime,dt=dt,ngear=ngear,block=NULL,mtype=mtype[p])
      dat=tmp$resid
      dat$t = (start.year-0.4*dt +  dat$t*dt)

      ncol = ncol(dat)
   #  print(head(dat))
      if (mtype[p] == "x")
      {
         gear.col = 6
         dd = c(2:3,(gear.col+1):ncol)
         y = matrix(0.0,nrow=ntime,ncol=3)
         dat[,dd] = exp(dat[,dd])
         y[,1] = dat$pop1
         y[,2] = dat$pop2
         y[,3] = dat$pop1 + dat$pop2
      } 
      else if (mtype[p] == "i")
      {
         gear.col = 4
         dd = c(2,(gear.col+1):ncol)
         dat[,dd] = exp(dat[,dd])
         y = matrix(0.0,nrow=ntime,ncol=1)
         y[,1] = dat$pop
      }
      else
         stop(paste("Unknown model type (",mtype[p],
                    ") passed to plot.biomass.array(...)",sep=""))

      dev.set(biomass.dev)
      plot.biomass(dat$t,y,sd=tmp$sdlogPop,K=dat$K, propL=dat$propL,
                   forcing=dat$forcing,yrange=yrange)
      legend(x="topleft",legend=legend,bty='n',cex=1.6)

      F.ndx=grep("F",names(dat)) 
      predC.ndx=grep("predC",names(dat)) 
      obsC.ndx=grep("obsC",names(dat)) 

      Fmort = rowSums(dat[,F.ndx])
      obsC = rowSums(dat[,obsC.ndx])
      predC = rowSums(dat[,predC.ndx])

      dev.set(prod.dev)
      plot.production(Fmort,obsC,predC,dat$t, r=tmp$r,K=tmp$K)
      legend(x="topleft",legend=legend,bty='n',cex=1.6)



   }

   dev.set(biomass.dev)
   din = par("din")
   save.png.plot("biomass-array",width=din[1],height=din[2])

   dev.set(prod.dev)
   din = par("din")
   save.png.plot("production-array",width=din[1],height=din[2])

#  par(old.par)
}

## diagnostic.array=function(file,ntime=61,dt=1,ngear=5,mtype=NULL,block=NULL)
## {
##    require(grid)
## #  require(gridBase)
## 
##    print(paste("Scanning file",file))
##    log = scan(file,what="character")
##    res = grep("Residuals:",log)
##    if (is.null(block))
##    {
##       block = length(res)
##    }
##    tmp=get.diagnostics(log,ntime=ntime,dt=dt,ngear=ngear,block=block,mtype=mtype)
##    print(names(tmp))
##    dat=tmp$resid
##    print(head(dat))
##    predC.ndx=grep("predC",names(dat)) 
##    obsC.ndx=grep("obsC",names(dat)) 
##    if (mtype == "x")
##    {
##       gear.col = 6
##       y = matrix(0.0,nrow=ntime,ncol=3)
##       y[,1] = dat$pop1
##       y[,2] = dat$pop2
##       y[,3] = dat$pop1 + dat$pop2
##    }
##    else if (mtype == "i")
##    {
##       gear.col = 4
##       y = matrix(0.0,nrow=ntime,ncol=1)
##       y[,1] = dat$pop
##    }
##    else
##       stop(paste("Unknown model type(",mtype,
##                  ") passed to diagnostic,array(...)",sep=""))
## 
##    width = 9.0
##    height = 6.5
## #  x11(width=width,height=height)
## #  par(mar=c(4,4,0,0)+0.1)
##    opar <- par(no.readonly=TRUE) # Saving graphical parameters
## #  la = layout(matrix(c(1,1,2,3),ncol=2,nrow=2,byrow=FALSE))
## #  layout.show(la)
## #  plot.catches(dat$t,dat[,c(gear.col+2*ngear+1):(gear.col+3*ngear)],
## #                     dat[,c(gear.col+ngear+1)  :(gear.col+2*ngear)],
## #                     tmp$sdlogYield)#,block) 
##    plot.new()
##    multitop.vp <- viewport(layout=grid.layout(2,2), width = unit(0.95, "npc"))
##    grid.show.layout(grid.layout(1,2)) 
##    pl1 <- viewport(layout.pos.col=1, layout.pos.row=2, name="A")
##    p21 <- viewport(layout.pos.col=2, layout.pos.row=1, name="B")
##    p22 <- viewport(layout.pos.col=2, layout.pos.row=2, name="C")
##    vpall <- vpTree(multitop.vp, vpList(pl1,p21,p22))
##    pushViewport(vpall)
## 
##    upViewport()
##    downViewport("A")
## # stackedplot(data=list(1:10,10:1,rep(10,10)),main="A")
## #  plot.biomass(dat$t,y,sd=tmp$sdlogPop,K=dat$K, forcing=dat$forcing)
## #  plot.catches(dat$t,dat[,c(gear.col+2*ngear+1):(gear.col+3*ngear)],
## #                     dat[,c(gear.col+ngear+1)  :(gear.col+2*ngear)],
## #                     tmp$sdlogYield)#,block) 
## 
##    upViewport()
##    downViewport("B")
## # stackedplot(data=list(10:1,rep(10,10),1:10),main="B")
## #  plot.biomass(dat$t,y,sd=tmp$sdlogPop,K=dat$K, forcing=dat$forcing)
##    
##    upViewport(2)
## 
##   
## 
## #  plot.biomass(dat$t,y,sd=tmp$sdlogPop,K=dat$K, forcing=dat$forcing)
## 
##  
##    par(opar) # Returning the graphical parameters saved earlier
## }
## 
## # opar <- par(no.readonly=TRUE) # Saving graphical parameters
## # plot.new() # Needed for par(new=TRUE) in stackedplot()
## # 
## # multitop.vp <- viewport(layout=grid.layout(1,2), width = unit(0.95, "npc"))
## # pl1 <- viewport(layout.pos.col=1, layout.pos.row=1, name="A")
## # pl2 <- viewport(layout.pos.col=2, layout.pos.row=1, name="B")
## # vpall <- vpTree(multitop.vp, vpList(pl1,pl2))
## # pushViewport(vpall)
## # 
## # upViewport()
## # downViewport("A")
## # stackedplot(data=list(1:10,10:1,rep(10,10)),main="A")
## # 
## # upViewport()
## # downViewport("B")
## # stackedplot(data=list(10:1,rep(10,10),1:10),main="B")
## # 
## # upViewport(2)
## # par(opar) # Returning the graphical parameters saved earlier
