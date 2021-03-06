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
#  print(paste("block number =",block.number))
   if (block.number > 1)
      mtext(text=paste("(",block.number,")",sep=""),side=1,line=line,
            at=c(x,0),cex=0.8,col=col)
}

plot.catches=function(t,obs,pred,sd=NULL,block=NULL)
{
   print(head(t))
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

plot.biomass=function(x,y,K=NULL,sd=NULL, sdI=NULL, block=NULL,propL=NULL,
                      forcing=NULL,yrange=NULL,B1=NULL,indexed=NULL)
{
   print("plot.biomass, diagnostics.R")
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
   print(paste("ncol",ncol(y)))

   if (!is.null(K))
   {
      lines(x,K,lwd=2,lty="dotdash",col="blue",xlim=xrange)
      text(x[1],K[1],"K ",adj=c(1.25,0.5),col="blue",cex=1.2)
   }

#  sdlogNN = sqrt(4.0*sd*sd) # this is probably not correct
#  print("N3 error")
#  if ( (ncol(y) > 2) && (!is.null(sd)) )
#     plot.error(x,y[,3],sdlogNN, bcol="blue",fcol="lightblue")
   print(paste("Index error sd =",sdI))
   if (!is.null(sdI))
   {
   #  print("Index error")
      plot.error(x,y[,1],sdI, bcol="purple4",fcol="purple")
   }

   print(paste("N error,sd =",sd))
   if (!is.null(sd))
   {
   #  print("N1 error")
      plot.error(x,y[,1],sd, bcol="blue",fcol="lightblue")
   }
   if ( (ncol(y) > 1) && (!is.null(sd)) )
   {
      print("N2 error")
      plot.error(x,y[,2],sd, bcol="blue",fcol="lightblue")
      lines(x,y[,1],col="blue",lwd=5)
   }
   lines(x,y[,1],col="blue",lwd=5)
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
#  if ( (!is.null(forcing)) && (max(forcing)> 1) )
   print(paste("p.b indexed",indexed))
   if ( (!is.na(indexed)) && (indexed) )
   {
      print("forcing")
      par("new"=TRUE)
   #  frange = yrange #c(0,1.2*max(forcing,na.rm=TRUE))
   #  ftic <- pretty(c(frange[1],frange[2]),5)
   #  plot(xrange,yrange,type='n',axes=FALSE,ann=FALSE,yaxs='i')
      lines(x,forcing,lwd=3,col="purple")
   #  plot(x,forcing,lwd=3,type='l',col="purple",
   #       ylim=frange,ann=FALSE,axes=FALSE,xlim=xrange)
      text(x[1],forcing[1],"I ",adj=c(1.5,0.5),col="purple",cex=1.2)
   #  axis(4,lwd=1,at=ftic[2:length(ftic)],col="purple",col.axis="purple")
   #  abline(v=par("usr")[2],lwd=2,col="purple")
   #  mtext("Abundance Index",side=4,col="purple",line=3)
    }
    if (!is.null(B1))
    {
       points(x[1],B1,pch=8,cex=1.2,col="blue")
    #  text(x[1],B1,adj=c(1.5,0.5),col="blue",cex=1.2)
       text(x[1],B1,parse(text=paste("B","[1]")),adj=c(1.25,0.5),col="blue",cex=1.2)
    }
    if (!is.null(block))
       show.block.number(block,x[1])
}


plot.production=function(Fmort, obsC, predC, t, r, K, block=NULL)
{
    print(paste("plot.prod",r,K))
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


read.rep.diagnostics=function(file="issams6.rep", ntime=61,dt=1,ngear=4,
                              block=NULL, mtype='i',klingon=FALSE)
{

   print(paste("Scanning file",file))
   rep = scan(file,what="character")
   print("head(rep):")
   print(head(rep))

   ests=get.diagnostics(rep,ntime=ntime,dt=dt,ngear=ngear,block=NULL,klingon=klingon)
   return(ests)
}

extract.value=function(text,text.list,sig.fig=3)
{
   name = paste("^",text,sep="")
#  print(paste(text,name))
   ndx = grep(name,text.list)[1]
#  print(paste("ndx",ndx))
#  print(paste("nchar",nchar(name)))
   colon = grep(":",name)
#  if (text=="K")
#  {
#     print(paste("ndx",ndx))
#     print(paste("colon",colon,length(colon),is.numeric(colon),as.numeric(colon),sep="|"))
#  }
   offset = 2
   if(length(colon)>0)
      offset = 1
#  print(paste("offset",offset))
#  print(text.list[ndx])
#  print(text.list[ndx+offset])

   x = as.numeric(text.list[ndx+offset])
   x = signif(x,sig.fig)

#  x = text.list[ndx+offset]

#  print(paste(text,x))
   return(x)
}


get.diagnostics=function(log,ntime=61,dt=1,ngear=4,block=NULL,mtype='i',klingon=FALSE)
{

#  print(paste("Scanning file",file))
#  log = scan(file,what="character")
## print(length(log))
   if (klingon)
      ngear = ngear + 1
   res = grep("Residuals:",log)
   nblock = length(res)
   print(paste(nblock,"residual blocks found"))

   if (is.null(block))
      block = nblock

   ests = list()
   ests$nll = extract.value("nll",log,sig=5)
   ests$nvar = extract.value("nvar",log)
   ests$Gmax = extract.value("Gmax",log)

   ests$logMSY = extract.value("logMSY",log)
   ests$logFmsy = extract.value("logFmsy",log)
   ests$logitQ = extract.value("logitQ:",log)
   ests$logsdlogProc = extract.value("logsdlogProc:",log)
   ests$logsdlogBProc = extract.value("logsdlogBProc:",log)
   ests$logsdlogN = extract.value("logsdlogN:",log)
   ests$logsdlogYield = extract.value("logsdlogYield:",log)

   ests$MSY = extract.value("MSY",log)
   ests$Fmsy = extract.value("Fmsy",log)
   ests$Q = extract.value("Q:",log)
   if (is.na(ests$logsdlogBProc))
      ests$logsdlogBProc = ests$logsdlogN
   if (is.na(ests$logsdlogProc))
      ests$logsdlogProc = ests$logsdlogBProc

   ests$sdlogBProc = signif(exp(ests$logsdlogBProc),3)
   ests$sdlogFmort = extract.value("sdlogFmort:",log)
   ests$sdlogYield = extract.value("sdlogYield:",log)
   ests$sdlogI     = extract.value("sdlogI:",log)

   ests$r = extract.value("r",log)
   ests$K = extract.value("K",log)

   ests$pcon = extract.value("pcon",log)
   ests$rprior = extract.value("r_prior",log)
   ests$sdrprior = extract.value("sdr_prior",log)
#  print(ests)

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

   # correct indexing for ngear
   Fndx=c(seq(5,(4+ngear)))
   print(Fndx)
   predCndx=c(seq((5+ngear),(4+2*ngear)))
   obsCndx=c(seq((5+2*ngear),(4+3*ngear)))
   print(predCndx)

   totalF = rowSums(exp(resid[,Fndx]))
   print(paste("ntime:",ntime))
   ests$lastF = signif(mean(totalF[seq(ntime-4,ntime)]),3)
   print(paste("Last 5 year total average F:",ests$lastF))
   totalC = rowSums(exp(resid[,obsCndx]))
   ests$lastC = signif(mean(totalC[seq(ntime-4,ntime)]),3)
   print(paste("Last 5 year average total catch:",ests$lastC))

   ests$averageC = signif(mean(totalC),3)
   print(paste("Average total catch:",ests$averageC))

#  print(paste("Block",block,":"))
#  print(ests)
#  print(head(resid))

   ests$resid=as.data.frame(resid)

   return(ests)
 
}

plot.diagnostic.array=function(file="issams.rep",ntime=61,dt=1,ngear=4)
{
   rep = read.rep.diagnostics(file=file)
   print(names(rep))
   dat = rep$resid
   print(head(dat))

   width = 9.0
   height = 9.0

   x11(width=width,height=height,title="Biomass",xpos=100)
   old.par = par(no.readonly = TRUE,"new"=FALSE) 
   par(mar=c(3,5,0,4)+0.1,"new"=FALSE)
   lm = layout(matrix(c(1,1,2,3),ncol=2,byrow=FALSE))

   layout.show(lm)

   gear.col = 4
   plot.catches(dat$t,dat[,c(gear.col+2*ngear+1):(gear.col+3*ngear)],
                      dat[,c(gear.col+ngear+1)  :(gear.col+2*ngear)],
                      rep$sdlogYield)#,block) 
   y = matrix(nrow=ntime,ncol=1)
   y[,1]=dat$pop
#  y[,2]=dat$forcing
   plot.biomass(dat$t,y,sd=rep$sdlogPop,block=NULL,K=dat$K,B1=B1,
                forcing=dat$forcing,indexed='i')
 


   par(old.par,"new"=FALSE) 
}



##plot.biomass.array=function(path.list=c("../run-issams-dev/issams-dev.rep",
##                                         "../run-issams/issams.rep",
##                                         "../run-xssams/1/xssams.rep"),
##plot.biomass.array=function(path.list=c("./run-issams/r2/Q0/issams.rep",
##                                         "./run-issams-dev/r2/Q0/issams-dev.rep",
##                                         "./run-issams/r2/Q1/issams.rep",
##                                         "./run-issams-dev/r2/Q1/issams-dev.rep",
##                                         "./run-xssams/1/xssams.rep"),
##                             ntime=61,dt=1,ngear=5,mtype=c("i","i","i","i","x"))
#plot.biomass.array = function(path.list=c("./run-issams/use_r_prior/r2/Q0/issams.rep",
#           "./run-issams/use_r_prior/r2/Q1/issams.rep",
#           "./run-issams-dev/use_r_prior/r2/Q0/issams-dev.rep",
#           "./run-issams-dev/use_r_prior/r2/Q1/issams-dev.rep"),
#           ntime=61,dt=1,ngear=5,mtype=c("i","i","i","i"))
#{
#   npath = length(path.list)
#   print(npath)
#   width = 15.0
#   height = 9.0
##  yrange=c(0,30000)
#   start.year = 1952
#
#   old.par = par(no.readonly = TRUE,"new"=FALSE) 
#
#   x11(width=width,height=height,title="Biomass",xpos=100)
#   par(mar=c(3,5,0,4)+0.1,"new"=FALSE)
#   lm = layout(matrix(c(1:npath),ncol=2,byrow=FALSE))
#   layout.show(lm)
#   biomass.dev = dev.cur()
#
#   x11(width=width,height=height,title="Production",xpos=200)
#   par(mar=c(4,4.5,0,0)+0.1,"new"=FALSE,las=1)
#   lm = layout(matrix(c(1:npath),ncol=2,byrow=FALSE))
#   layout.show(lm)
#   prod.dev = dev.cur()
#
#   for (p in 1:npath)
#   {
#      path=path.list[p]
#      print(path)
#      pv = unlist(strsplit(path,"[\\./]"))
#      model = pv[length(pv)-1]
#      print(model)
#   #  if (model == "issams-dev")
#   #     #     parse(text=paste("N","[1]","+","N","[2]"))
#   #     legend = parse(text=paste("B","[1]","~d","~~no~index"))
#   #  else if (model == "issams")
#   #     legend = parse(text=paste("MSY~F","[msy]","~~indexed"))
#   #  else if (model == "xssams")
#   #     legend = parse(text=paste("MSY~F","[msy]","~~immigration"))
#   #  else
#   #     legend = model
#      if (p == 1)
#         legend = parse(text=paste("A.~MSY~F","[msy]","~~no~index"))
#      else if (p == 2)
#         legend = parse(text=paste("C.~MSY~F","[msy]","~~indexed"))
#      else if (p == 3)
#         legend = parse(text=paste("B.~B","[1]","~d","~~no~index"))
#      else if (p == 4)
#         legend = parse(text=paste("D.~B","[1]","~d","~~indexed"))
#      else
#         legend = model
#
#      log = scan(path,what="character")
#      res = grep("Residuals:",log)
#      block = length(res)
#      tmp=get.diagnostics(log,ntime=ntime,dt=dt,ngear=ngear,block=NULL,mtype=mtype[p])
#      print(paste("tmp$indexed",tmp$indexed))
#      dat=tmp$resid
#      dat$t = (start.year-0.4*dt +  dat$t*dt)
#
#      ncol = ncol(dat)
#   #  print(head(dat))
#      if (mtype[p] == "x")
#      {
#         gear.col = 6
#         dd = c(2:3,(gear.col+1):ncol)
#         y = matrix(0.0,nrow=ntime,ncol=3)
#         dat[,dd] = exp(dat[,dd])
#         y[,1] = dat$pop1
#         y[,2] = dat$pop2
#         y[,3] = dat$pop1 + dat$pop2
#      } 
#      else if (mtype[p] == "i")
#      {
#         gear.col = 4
#         dd = c(2,(gear.col+1):ncol)
#         dat[,dd] = exp(dat[,dd])
#         y = matrix(0.0,nrow=ntime,ncol=1)
#         y[,1] = dat$pop
#      }
#      else
#         stop(paste("Unknown model type (",mtype[p],
#                    ") passed to plot.biomass.array(...)",sep=""))
#
#      dev.set(biomass.dev)
#   #  print(paste("tmp$indexed",tmp$indexed))
#      yrange = c(0,1.2*max(y,na.rm=TRUE,dat$K))
#      plot.biomass(dat$t,y,sd=tmp$sdlogPop,K=dat$K, propL=dat$propL,
#                   forcing=dat$forcing,B1=tmp$B1,yrange=yrange,
#                   indexed=tmp$indexed)
#      legend(x="topleft",legend=legend,bty='n',cex=1.6)
#
#      F.ndx=grep("F",names(dat)) 
#      predC.ndx=grep("predC",names(dat))
#      obsC.ndx=grep("obsC",names(dat)) 
#
#      Fmort = rowSums(dat[,F.ndx])
#      obsC = rowSums(dat[,obsC.ndx])
#      predC = rowSums(dat[,predC.ndx])
#
#      dev.set(prod.dev)
#      plot.production(Fmort,obsC,predC,dat$t, r=tmp$r,K=tmp$K)
#      legend(x="topleft",legend=legend,bty='n',cex=1.6)
#
#
#
#   }
#
#   dev.set(biomass.dev)
#   din = par("din")
#   save.png.plot("biomass-array",width=din[1],height=din[2])
#
#   dev.set(prod.dev)
#   din = par("din")
#   save.png.plot("production-array",width=din[1],height=din[2])
#
#   par(old.par,"new"=FALSE) 
#}
#
#
#
#make.fit.table = function(paths=c("./run-issams/use_r_prior/r2/Q0/issams",
#           "./run-issams-dev/use_r_prior/r2/Q0/issams-dev",
#           "./run-issams/use_r_prior/r2/Q1/issams",
#           "./run-issams-dev/use_r_prior/r2/Q1/issams-dev"))
#{
##  anames=c("alogB1", "alogdB1K", "aB1", "adB1K", "aMSY", "aFmsy", "alogr",
##           "ar", "aK", "asdlogProc", "asdlogYield", "aQ", "alogQ", "apcon")
#   anames=c("aB1", "adB1K", "aMSY", "aFmsy", "ar", "aK", 
#             "asdlogProc", "asdlogYield", "aQ","aT12","aT21")
#   tnames=c("$B_1$", "$d$", "$\\MSY$", "$\\Fmsy$", "$r$", "$K$", 
#             "$\\sigma_P$", "$\\sigma_Y$", "$Q$","$T_{12}$","$T_{21}$")
#   nnames = length(anames)
##  indexing=c("Q0","Q1")
##  models=c("issams","issams-dev")
##  nmod = length(indexing)*length(models)
##  paths=c("./run-issams/r2/Q0/issams",
##          "./run-issams-dev/r2/Q0/issams-dev",
##          "./run-issams/r2/Q1/issams",
##          "./run-issams-dev/r2/Q1/issams-dev",
##          "./run-xssams/r2/xssams")
#
## AIC = -2 ln L +2p and BIC = -2 ln L + p ln(n) 
## p = number of parameters
## n = number of observations
#
#
#   varnames=c("$n$","$-\\log L$","$|G|_{max}$","AIC",tnames)
#   fits = matrix(ncol=length(paths),nrow=length(varnames))
#   rownames(fits) = varnames
#   
#   mm = 0
##  for (Q in indexing)
##  {
##     for (m in modelsr
#      for (p in paths)
#      {
#         mm = mm + 1
#      #  path=paste(root,"/run-",m,"/",region,"/",Q,sep="")
#      #  print(paste(mm,path))
#      #  par.file = paste(path,"/",m,".par",sep="")
#      #  std.file = paste(path,"/",m,".std",sep="")
#         par.file = paste(p,".par",sep="")
#         std.file = paste(p,".std",sep="")
# 
#         print(par.file)
#         par = scan(file=par.file,what="character")
#         npar = as.numeric(par[6])
#         nll  = -1*as.numeric(par[11])
#         gmax = as.numeric(par[16])
#         AIC = 2*(npar - as.numeric(par[11]))
#         print(paste(npar,nll,gmax))
#      #  fits[1,mm] = m
#      #  fits[2,mm] = Q
#         fits[1,mm] = npar
#         fits[2,mm] = nll
#         fits[3,mm] = gmax
#         fits[4,mm] = AIC
#         r = 4
#
##  print(std.file)
#         print(std.file)
#         std=try(scan(file=std.file,what="character"),TRUE)
#         print(class(std))
#         if (class(std) == "try-error")
#         {
#            print(paste(std.file,"missing"))
#         }
#         else
#         {
#            for (i in 1:nnames)
#            {
#               w = which(std == anames[i])
#               print(paste(i,anames[i],w,length(w),std[w+1],std[w+2]))
#               if ( (length(w)>0) && (as.numeric(std[w+2] > 1.0e-3)) )
#                  fits[r+i,mm] = as.numeric(std[w+1])
#            }  
#         }
#      #  print(head(fits))
#      }
##  }
#
#   write.table(fits, file = "fit_table.csv", append = FALSE, quote = FALSE, sep = ",",
#               eol = "\n", na = "NA", dec = ".", col.names = FALSE,
#               row.names = varnames) #, qmethod = c("escape", "double"),
#             # fileEncoding = "")
#
#
#
#   return(fits)
#}


print.error=function(y=0.486,sd=0.8,mult=2)
{
   #  sdyu = exp(log(y)+mult*sd)
   #  sdyl = exp(log(y)-mult*sd)
   yu = exp(log(y)+mult*sd)
   yl = exp(log(y)-mult*sd)

   print(paste(y,sd,yl,yu))
#  xx = seq(y/5,y*5,0.01)
   xx = seq(0.01,2.5,0.01)
#  print(xx)
   yy = dnorm(log(xx),mean=log(y),sd=sd)
   plot(xx,yy,type='l')
   abline(v=y,col="blue")
   abline(v=yl,col="blue",lty="dotdash")
   abline(v=yu,col="blue",lty="dotdash")
}

plot.obs.pred.catch=function(file="issams.rep", ntime=61,dt=1,ngear=4,
                             block=NULL, mtype='i')
{
   rep = read.rep.diagnostics(file=file, ntime=ntime, dt=dt,
                              ngear=ngear, block=blocak, mtype=mtype)

   pred = exp(rep$resid[9:12])
   obs = exp(rep$resid[13:16])
   cmax = max(obs,pred)
   print(cmax)

   plot(c(0,cmax),c(0,cmax),xlab="Observed (mt)",ylab="Predicted (mt)",type='l',
        axes=FALSE,main="Observed and Predicted Catch")
   axis(1,lwd=0)
   abline(h=par("usr")[3],lwd=3)
   axis(2,lwd=0)
   abline(v=par("usr")[1],lwd=3)

   for (g in 1:ngear)
   {
      points(obs[,g],pred[,g])#,pch=16,col=g+1)
   }
      

   save.png.plot("OPcatch",width=6.5,height=6.5)
}

#plot.multi.obs.pred.catch=function(path.list=c("r2","r4","r0-sdrprior"), 
plot.multi.obs.pred.catch=function(path.list=c("r2","r2-sdrprior"), 
                                   ntime=61,dt=1,ngear=4, block=NULL, mtype='i')

{
   cmax = 1000
   plot(c(0,cmax),c(0,cmax),xlab="Observed (mt)",ylab="Predicted (mt)",type='l',
        axes=FALSE,main="Observed and Predicted Catch")
   axis(1,lwd=0)
   abline(h=par("usr")[3],lwd=3)
   axis(2,lwd=0)
   abline(v=par("usr")[1],lwd=3)
   pcol = 1
   ppch = 1
   clist = c("blue","red") 
   for(p in path.list)
   {
      if (length(grep("2",p)) > 0)
         ppch = 1 # o
      else
         ppch = 3 # +

      if (nchar(p) > 2)
         pcol = 2 # red
      else
         pcol = 1 # blue
         
      file = paste(p,"/issams6.rep",sep="")
      rep = read.rep.diagnostics(file=file, ntime=ntime, dt=dt,
                              ngear=ngear, block=blocak, mtype=mtype)
      # hard coded indices; only valid for ngear=4
      pred = exp(rep$resid[9:12])
      obs = exp(rep$resid[13:16])
      print(max(obs,pred))

      for (g in 1:ngear)
      {
      #  points(obs[,g],pred[,g],pch=ppch)
         points(obs[,g],pred[,g],col = clist[pcol],pch=ppch) #,col=g+1)
      }
   #  ppch = 3
      
   }
   save.png.plot("multiOPcatch",width=6.5,height=6.5)
}
