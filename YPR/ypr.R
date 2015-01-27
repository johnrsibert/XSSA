# the functions in this file use functions in utilities.R and ReadADMB.R
# source("../scripts/utils.R")
# source("../scripts/ReadADMB.R")

#also:
require("R4MFCL")


#from plot-12.par.rep

# Beverton-Holt stock-recruitment relationship report
# alpha = 1.418e+008  beta =  1.370e+005  steepness = 8.000e-001

#from 12.par
# The von Bertalanffy parameters
# 25.504923890770  20.000000000000   40.000000000000 
# 153.437994496285  140.000000000000   200.000000000000 
# 0.137527215929   0   0.300000000000 

#from yft.ini
# The von Bertalanffy parameters
# Initial  lower bound  upper bound
# ML1
# 25 20 40
# ML2
# 150 140 200
# K (per year)
# 0.15 0 0.3
# Length-weight parameters
# 2.512e-05 2.9396
mfcl.wa = 2.512e-05
mfcl.wb = 2.9396

# from Dalzell
# The current minimum size for YFT for commercial sale in Hawaii is 3 lb.
# Is it possible that the Y/R analysis could look at a range of potential minimum sizes.
# I believe Dave Itano has suggested 10 lb, but maybe we would want to also check out 
# what would be the impact of 15 lb and 20 lb on the Y/R.

kgperlb = 0.453592
mfcl.rep.file = paste(path.expand("~"),
                     "/Projects/xssa/mfcl/yft/2014/basecase/plot-12.Rdata",
                     sep="")


ypr.comp=function(WatAge,MatAge,FatAge)
{
   nages = length(WatAge)
   Y = vector(length=nages)
   prevW = WatAge[1]
   prevFF = FatAge[1]
   prevN = 1.0
   prevYY = prevFF*prevN*prevW
   Y[1] = prevYY
   for (a in 2:nages)
   {
      W = WatAge[a]
      FF = FatAge[a]
      Z = FF+MatAge[a]
      N = prevN*exp(-Z)	
      YY = FF*N*W

      prevW = W
      prevFF = FF
      prevN = N
      prevYY = YY
      Y[a] = YY
   }
   YPR = sum(Y)
   return(YPR)
}


mfcl.ypr<-function(epoch=NULL,region=2,astar=NULL)
{

   rep.file = mfcl.rep.file
   print(rep.file)
   load(rep.file)

   nTime = rep$nTime
   nAges = rep$nAges
   nReg = rep$nReg
   yrs = rep$yrs
   MatAge = rep$MatAge
   MeanWatAge = rep$mean.WatAge
   MeanLatAge = rep$mean.LatAge
   FatAgeReg = rep$FatYrAgeReg
   if (is.null(astar))
      astar=nAges
   print(paste("astar =",astar))

   if (is.null(epoch))
      epoch = (nTime-19):nTime
   print(paste("Averaging over",length(epoch),"quarter epoch:"))
   print(epoch)

   AveFatAge = matrix(0,nrow=nReg,ncol=nAges)

   for (i in 1:nReg)
   {
      den = 0
      for (j in epoch)
      {
         den = den + 1
         for (k in 1:nAges)
         {
            AveFatAge[i,k] = AveFatAge[i,k] + FatAgeReg[j,k,i]
         }
      } 
      for (k in 1:nAges)
      {
         AveFatAge[i,k] = AveFatAge[i,k]/den
      }
   }

   FR = AveFatAge[region,]

   Y = vector(length=nAges)
   Fmult = seq(0.1,10,.1)
   nf = length(Fmult)
   YPR.f = vector(length=nf)
#  m = 10
   for (m in 1:nf)
   {
      FF = AveFatAge[region,]*Fmult[m]
      YPR.f[m] = ypr.comp(MeanWatAge,MatAge,FF)
   }

   YPR.a = vector(length=nAges)
#  m = 10
   maxa = 0
   maxYPR = 0
   for (m in 1:nAges)
   {
      FF = AveFatAge[region,]
      for (a in 1:m)
         FF[a] = 0.0;
      YPR.a[m] = ypr.comp(MeanWatAge,MatAge,FF)
      if (YPR.a[m] > maxYPR)
      {
         maxYPR = YPR.a[m]
         maxa = m;
      }
   }
   print(paste(maxa,MeanWatAge[maxa],MeanWatAge[maxa]/kgperlb))


   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,4,0,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   lm <- layout(matrix(c(1:2),nrow=2,byrow=TRUE))
   layout.show(lm)

   nice.ts.plot(Fmult,YPR.f,
         xlab="F multiplier",ylab="Yield per Recruit (kg)",lwd=7)
   abline(v=1.0,col="red",lty="dotdash")
   label.panel("A")

   nice.ts.plot(MeanWatAge,YPR.a,
         xlab="Weight at First Capture (kg)",ylab="Yield per Recruit (kg)",lwd=7)
   points(c(3*kgperlb,10*kgperlb,15*kgperlb,20*kgperlb),c(0,0,0,0),col="red",pch='|',cex=2)
   abline(v=MeanWatAge[maxa],col="red",lty="dotdash")
   label.panel("B")

   save.png.plot(paste("YPR_MFCL_R",region,sep=""),width=width,height=height)
   par(old.par)
   return(as.data.frame(cbind(YPR.f,Fmult)))
}


csm.ypr=function()
{
   rep = csm.read.rep()
   print(names(rep))
   fit = read.fit("csm")

   nmort = fit$npar/4
   qq = nmort*2
   qF = (qq+1):(qq+nmort)
   qM = (qq+nmort+1):(qq+2*nmort)

   nAges = nmort
   astar = nAges

   Fmult = seq(0.1,10,.1)
   nf = length(Fmult)
   YPR.f = vector(length=nf)
#  m = 10
   for (m in 1:nf)
   {
      FF = fit$est[qF]*Fmult[m]
      YPR.f[m] = ypr.comp(rep$mort.mat$wt,fit$est[qM],FF)
   }

   YPR.a = vector(length=nAges)
#  m = 10
   maxa = 0
   maxYPR = 0
   for (m in 1:nAges)
   {
      FF= fit$est[qF]
      for (a in 1:m)
         FF[a] = 0.0;
      YPR.a[m] = ypr.comp(rep$mort.mat$wt,fit$est[qM],FF)
      if (YPR.a[m] > maxYPR)
      {
         maxYPR = YPR.a[m]
         maxa = m;
      }
   }
   print(paste(maxa,rep$mort.mat$wt[maxa],rep$mort.mat$wt[maxa]/kgperlb))

   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,4,0,0)+0.1)
   lm <- layout(matrix(c(1:2),nrow=2,byrow=TRUE))
   layout.show(lm)

   nice.ts.plot(Fmult,YPR.f,
         xlab="F multiplier",ylab="Yield per Recruit (kg)",lwd=7)
   abline(v=1.0,col="red",lty="dotdash")
   label.panel("A")

   nice.ts.plot(rep$mort.mat$wt,YPR.a,
         xlab="Weight at First Capture (kg)",ylab="Yield per Recruit (kg)",lwd=7)
   points(c(3*kgperlb,10*kgperlb,15*kgperlb,20*kgperlb),c(0,0,0,0),col="red",pch='|',cex=2)
   label.panel("B")
   abline(v=rep$mort.mat$wt[maxa],col="red",lty="dotdash")

   save.png.plot("YPR_csm",width=width,height=height)
   par(old.par)

}
csm.read.rep=function(path=".")
{
   get.field<-function(sca)
   {
     p = .csm.field.counter
     field = sca[p]
     p = p+1;
     .csm.field.counter<<-p
     return(field)
   }
   get.numeric.field<-function(sca,nfield)
   {
      ret<-as.numeric(get.field(sca))
      return(ret)
   }

   #######################################

   fp = paste(path,"csm.rep",sep='/')
   print(paste("Scanning ",fp))
   sca = scan(file=fp,what="raw")
   print(paste("Scan read",length(sca),"fields"),quote=F)
   .csm.field.counter <<- 1
   ret = list()

   f = which(sca == "#ngroup_MF")
   .csm.field.counter <<- f+1
   ng = get.numeric.field(sca)
   ncol = 6
   .csm.field.counter = f + 8
   qM = vector(length=ng)
   qF = vector(length=ng)
   wt = vector(length=ng)
   len = vector(length=ng)
   for (g in 1:ng)
   {
      n = get.numeric.field(sca)
      len[n] = get.numeric.field(sca)
      wt[n] = get.numeric.field(sca)
      nobs = get.numeric.field(sca)
      qM[n] = get.numeric.field(sca)
      qF[n] = get.numeric.field(sca)
   }
   print(nobs)
   ret$mort.mat = as.data.frame(cbind(len,wt,qM,qF))


   f = which(sca == "days")
#  print(sca[(f+3):(f+5)])
   .csm.field.counter = f+3
   ntime = 101
   days=vector(length=ntime)
   obs=vector(length=ntime)
   pred=vector(length=ntime)
   for (t in 1:ntime)
   {
      days[t] = get.numeric.field(sca)
      obs[t] = get.numeric.field(sca)
      pred[t] = get.numeric.field(sca)
   }
   ret$attrition=as.data.frame(cbind(days,obs,pred))

  return(ret)
}

plot.mfcl.mortality=function(epoch=NULL)
{
   print(paste("Loading ",mfcl.rep.file))
   load(mfcl.rep.file)

   nTime = rep$nTime
   nAges = rep$nAges
   nReg = rep$nReg
   yrs = rep$yrs
   MatAge = rep$MatAge
   MeanWatAge = rep$mean.WatAge
   MeanLatAge = rep$mean.LatAge
   FatAgeReg = rep$FatYrAgeReg

   # compute recet average fishing mortality for region
   if (is.null(epoch))
      epoch = (nTime-19):nTime
   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,4.5,1,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   lm <- layout(matrix(c(1:3),nrow=3,byrow=TRUE))
   layout.show(lm)

   nice.ts.plot(MeanWatAge,MatAge,
        xlab="Weight (kg)",lwd=7,
        ylab=substitute(paste("Natural Mortality ",(q^{-1}))))
#  legend("topleft",bty='n',cex=1.6,legend="M")
   label.panel("M")

   for (r in c(2,4))
   {
      den = 0
      AveFatAge = vector(length=nAges)
      for (j in epoch)
      {
         den = den + 1
         for (k in 1:nAges)
         {
            AveFatAge[k] = AveFatAge[k] + FatAgeReg[j,k,r]
         }
      } 
      for (k in 1:nAges)
      {
         AveFatAge[k] = AveFatAge[k]/den
      }
      
      nice.ts.plot(MeanWatAge,AveFatAge,
           ylab=substitute(paste("Fishing Mortality ",(q^{-1}))),
           xlab="Weight (kg)",lwd=7)
      points(c(3*kgperlb,10*kgperlb,15*kgperlb,20*kgperlb),c(0,0,0,0),col="red",pch='|',cex=2)
  #   legend("topleft",bty='n',cex=1.6,legend=paste("F(",r,")",sep=""))
      label.panel(paste("F(",r,")",sep=""))
   }

   save.png.plot("MFCL_MF",width=width,height=height)
   par(old.par)
}


plot.csm.mortality=function(epoch=NULL)
{
   print(paste("Loading ",mfcl.rep.file))
   load(mfcl.rep.file)

   nTime = rep$nTime
   nAges = rep$nAges
   nReg = rep$nReg
   yrs = rep$yrs
   MatAge = rep$MatAge
   MeanWatAge = rep$mean.WatAge
   MeanLatAge = rep$mean.LatAge
   FatAgeReg = rep$FatYrAgeReg

   # compute recet average fishing mortality for region
   if (is.null(epoch))
      epoch = (nTime-19):nTime

   AveFatAge = matrix(0,ncol=2,nrow = nAges)
   region.vec = as.vector(c(2,4))
   for (r in 1:2)
   {
      den = 0
      region = region.vec[r]
      print(paste("region",region))
      for (j in epoch)
      {  
         den = den + 1
         for (k in 1:nAges)
         {
            AveFatAge[k,r] = AveFatAge[k,r] + FatAgeReg[j,k,region]
         }
      }
      for (k in 1:nAges)
      {
         AveFatAge[k,r] = AveFatAge[k,r]/den
      }
   } 
   print(head(AveFatAge))

   csm.rep = csm.read.rep()
#  print(csm.rep$mort.mat$wt)
#  print(csm.rep$mort.mat$qF)
   csm = read.fit("csm")
   nmort = csm$npar/4
   qq = nmort*2
   qF = (qq+1):(qq+nmort)
   print(csm$names[qF])
   qM = (qq+nmort+1):(qq+2*nmort)
   print(csm$names[qM])

   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,4,0,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   lm <- layout(matrix(c(1:2),nrow=2,byrow=TRUE))
   layout.show(lm)

   maxM = 1.2*max(MatAge,csm.rep$mort.mat$qM)
   print(maxM)
   nice.ts.plot(MeanWatAge,MatAge,ylim=c(0,maxM),bcol="black",fcol="lightgray",
                xlab="Weight (kg)", 
                ylab=substitute(paste("Natural Mortality ",(q^{-1}))))
   double.lines(csm.rep$mort.mat$wt,csm$est[qM],bcol="blue",fcol="lightblue",lwd=7)
   sd.bars(csm.rep$mort.mat$wt,csm$est[qM],2.0*csm$std[qM],lwd=5,col="blue")
   label.panel("A")

   maxF = 1.2*max(AveFatAge,csm$est[qF])
   print(maxF)
   nice.ts.plot(MeanWatAge,AveFatAge,bcol="black",fcol="lightgray",
                legend=c(" R2"," R4"), ylim = c(0,maxF),
                xlab="Weight (kg)", 
                ylab=substitute(paste("Fishing Mortality ",(q^{-1}))))
   double.lines(csm.rep$mort.mat$wt,csm$est[qF],bcol="blue",fcol="lightblue",lwd=7)
   sd.bars(csm.rep$mort.mat$wt,csm$est[qF],2.0*csm$std[qF],lwd=5,col="blue")
   label.panel("B")


   save.png.plot("csm_MF",width=width,height=height)
   par(old.par)
   
}

LL.nonLL.catch=function()
{
   dat = read.table("../run/five_gears.dat")
   nTime = ncol(dat)
   ngear = nrow(dat)
   epoch = (nTime-19):nTime
#  print(t(dat[,epoch]))
   sum = vector(length=2,mode="numeric")
   den = vector(length=2,mode="numeric")
   for (g in 1:ngear)
   {
      if (g != 3)
         n = 1
      else
         n = 2
      for (j in epoch)
      {
      #  print(paste(g,n,j,dat[g,j]))
         sum[n] = sum[n] + dat[g,j]
         den[n] = den[n] + 1
      }
   }
   ave = sum/5
#  print(sum)
#  print(den)
#  print(ave)
   print(paste("Average annual non-lonline catch:",ave[1]))
   print(paste("Average annual     lonline catch:",ave[2]))
   return(ave)
}

bogus.ypr<-function(bogus.mode=4,epoch=NULL,region=2,astar=NULL)
{

   rep.file = mfcl.rep.file
   print(rep.file)
   load(rep.file)

   nTime = rep$nTime
   nAges = rep$nAges
   nReg = rep$nReg
   yrs = rep$yrs
   MatAge = rep$MatAge
   MeanWatAge = rep$mean.WatAge
   MeanLatAge = rep$mean.LatAge
   FatAgeReg = rep$FatYrAgeReg
   if (is.null(astar))
      astar=nAges
   print(paste("astar =",astar))

   if (is.null(epoch))
      epoch = (nTime-19):nTime
   print(paste("Averaging over",length(epoch),"quarter epoch:"))
   print(epoch)

   AveFatAge = matrix(0,nrow=nReg,ncol=nAges)

   for (i in 1:nReg)
   {
      den = 0
      for (j in epoch)
      {
         den = den + 1
         for (k in 1:nAges)
         {
            AveFatAge[i,k] = AveFatAge[i,k] + FatAgeReg[j,k,i]
         }
      } 
      for (k in 1:nAges)
      {
         AveFatAge[i,k] = AveFatAge[i,k]/den
      }
   }

   FR = AveFatAge[region,]
   maxFR = max(FR)
   maxF4a = 0
   maxF4 = 0
   for (a in 1:nAges)
   {
      if (AveFatAge[4,a] > maxF4)
      {
         maxF4a = a
         maxF4  = AveFatAge[4,a] 
      }
   }

   average.catch = LL.nonLL.catch()
   catch.ratio = average.catch[1]/average.catch[2]

   fB = dlnorm(1:nAges,log(bogus.mode),log(2.0))
   print(paste(sum(fB),sum(FR)))
   sfB = fB/max(fB)*maxFR*catch.ratio
   sFR = FR
   print(paste(sum(sfB),sum(sFR)))
 
   bogusF = sFR + sfB


# do YPR
   
   Y = vector(length=nAges)
   Fmult = seq(0.1,10,.1)
   nf = length(Fmult)
   YPR.f = vector(length=nf)
#  m = 10
   for (m in 1:nf)
   {
      FF = bogusF*Fmult[m]
      YPR.f[m] = ypr.comp(MeanWatAge,MatAge,FF)
   }

   YPR.a = vector(length=nAges)
#  m = 10
   maxa = 0
   maxYPR = 0
   for (m in 1:nAges)
   {
      FF = bogusF
      for (a in 1:m)
         FF[a] = 0.0;
      YPR.a[m] = ypr.comp(MeanWatAge,MatAge,FF)
      if (YPR.a[m] > maxYPR)
      {
         maxYPR = YPR.a[m]
         maxa = m;
      }
   }
   print(paste(maxa,MeanWatAge[maxa],MeanWatAge[maxa]/kgperlb))


   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,4,0,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   lm <- layout(matrix(c(1:3),nrow=3,byrow=TRUE))
   layout.show(lm)

   nice.ts.plot(MeanWatAge,bogusF,
           ylab=substitute(paste("Fishing Mortality ",(q^{-1}))),
           xlab="Weight (kg)",lwd=7)
   lines(MeanWatAge,sFR,col="blue",lwd=3,lty="dotted")
   lines(MeanWatAge,sfB,col="red",lwd=3,lty="dotted")
   points(c(3*kgperlb,10*kgperlb,15*kgperlb,20*kgperlb),c(0,0,0,0),col="red",pch='|',cex=2)
   abline(v=MeanWatAge[bogus.mode],col="red",lty="dotdash")
   label.panel("A")

   nice.ts.plot(Fmult,YPR.f,
         xlab="F multiplier",ylab="Yield per Recruit (kg)",lwd=7)
   abline(v=1.0,col="red",lty="dotdash")
   label.panel("B")

   nice.ts.plot(MeanWatAge,YPR.a,
         xlab="Weight at First Capture (kg)",ylab="Yield per Recruit (kg)",lwd=7)
   points(c(3*kgperlb,10*kgperlb,15*kgperlb,20*kgperlb),c(0,0,0,0),col="red",pch='|',cex=2)
   abline(v=MeanWatAge[maxa],col="red",lty="dotdash")
   label.panel("C")

   save.png.plot(paste("YPR_BOGUS",bogus.mode,sep=""),width=width,height=height)
   par(old.par)
   return(as.data.frame(cbind(YPR.f,Fmult)))
}

   
do.it.all=function()
{
   graphics.off()
   plot.mfcl.mortality()->junk
   mfcl.ypr(region=2)->junk
   mfcl.ypr(region=4)->junk
   plot.csm.mortality()
   bogus.ypr()->junk
   csm.ypr()
   LL.nonLL.catch()
}

