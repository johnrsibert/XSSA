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

mfcl.ypr<-function(epoch=NULL,region=2,astar=NULL)
{

   rep.file = paste(path.expand("~"),
                     "/Projects/xssa/mfcl/yft/2014/basecase/plot-12.Rdata",
                     sep="")
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
      Y = vector(length=nAges)
      prevW = MeanWatAge[1]
      prevFF = AveFatAge[region,1]*Fmult[m]
      prevN = 1.0
      prevYY = prevFF*prevN*prevW
      Y[1] = prevYY
      for (a in 2:nAges)
      {
         W = MeanWatAge[a]
         if (a <= astar)
            FF = AveFatAge[region,a]*Fmult[m]
         else
            FF = AveFatAge[region,a]
         Z = FF+MatAge[a]
         N = prevN*exp(-Z)	
         YY = FF*N*W
   
         prevW = W
         prevFF = FF
         prevN = N
         prevYY = YY
         Y[a] = YY
      }
      YPR.f[m] = sum(Y)
   }

   YPR.a = vector(length=nAges)
#  m = 10
   maxa = 0
   maxYPR = 0
   for (m in 1:nAges)
   {
      Y = vector(length=nAges)
      prevW = MeanWatAge[1]
      prevFF = 0.0
      prevN = 1.0
      prevYY = prevFF*prevN*prevW
      Y[1] = prevYY
      for (a in 2:nAges)
      {
         W = MeanWatAge[a]
         if (a <= m)
            FF = 0.0
         else
            FF = AveFatAge[region,a]
         Z = FF+MatAge[a]
         N = prevN*exp(-Z)	
         YY = FF*N*W
   
         prevW = W
         prevFF = FF
         prevN = N
         prevYY = YY
         Y[a] = YY
      }
      YPR.a[m] = sum(Y)
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
#  nice.ts.plot(MeanWatAge,AveFatAge[region,],
#        ylab="Fishing Mortality",xlab="Weight (kg)")
#  points(c(3*kgperlb,10*kgperlb,15*kgperlb,20*kgperlb),c(0,0,0,0),col="red",pch='|',cex=2)
#  if(astar < nAges)
#     abline(v=MeanWatAge[astar],col="green",lty="dotdash")

   nice.ts.plot(Fmult,YPR.f,
         xlab="F multiplier",ylab="Yield per Recruit (kg)")
   abline(v=1.0,col="red",lty="dotdash")
   legend("topleft",bty='n',cex=1.6,legend="A")
#  title(main=paste("MFCL Region",region),line=-2,cex.main=1.6)

   nice.ts.plot(MeanWatAge,YPR.a,
         xlab="Weight at First Capture",ylab="Yield per Recruit (kg)")
   points(c(3*kgperlb,10*kgperlb,15*kgperlb,20*kgperlb),c(0,0,0,0),col="red",pch='|',cex=2)
   abline(v=MeanWatAge[maxa],col="red",lty="dotdash")
   legend("topleft",bty='n',cex=1.6,legend="B")

   save.png.plot(paste("YPR_MFCL_R",region,sep=""),width=width,height=height)
   par(old.par)
   return(as.data.frame(cbind(YPR.f,Fmult)))
}

mfcl.plot.yield<-function(epoch=NULL,region=2)
{
   rep.file = paste(path.expand("~"),
                     "/Projects/xssa/mfcl/yft/2014/basecase/plot-12.Rdata",
                     sep="")
   print(rep.file)

   load(rep.file)
#  print(names(rep))

   nTime = rep$nTime
   nAges = rep$nAges
   nReg = rep$nReg
   yrs = rep$yrs
   MatAge = rep$MatAge
   MeanWatAge = rep$mean.WatAge
   MeanLatAge = rep$mean.LatAge
   plot(1:nAges,MatAge,type='l')
#  print(paste(nAges,length(MatAge)))
   FatAgeReg = rep$FatYrAgeReg
#  print(dim(FatAgeReg))
   if (is.null(epoch))
      epoch = (nTime-19):nTime

   print(paste("Averaging over",length(epoch),"quarter epoch:"))
   print(epoch)

   AveFatAge = matrix(0,nrow=nReg,ncol=nAges)

   for (i in 1:nReg)
   {
#     print(paste("Region",i))
#     print(AveFatAge[i,])

      den = 0
   #  for (j in 1:20)
   #  for (j in (nTime-19):nTime)
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
#     print(paste("denominator =",den))
#     print(AveFatAge[i,])
   }

   x11()
   plot(1:nAges,AveFatAge[region,],type='l')
#  lines(1:nAges,AveFatAge[2,],col="blue")


   N = 1.0
   Y = 0.0
   nstep = 100
   dt = 1.0/nstep
   NatAge = vector(mode="numeric",length=nAges)
   YatAge = vector(mode="numeric",length=nAges)
   AveFatAge = AveFatAge
   for (k in 1:nAges)
   {
      Z = AveFatAge[region,k]+MatAge[k]
      Y = 0.0
      for(it in 1:nstep)
      {
	 N = N*exp(-Z*dt)
         Y = Y + AveFatAge[region,k]*N*MeanWatAge[k] * dt
      }
      NatAge[k] = N
      YatAge[k] = Y
   }
   print(sum(YatAge))
#  YatAge = YatAge/sum(YatAge)
 
   x11()
   plot(1:nAges,YatAge,type='l')
#  plot(1:nAges,YatAge,type='l',axes=FALSE)
#  axis(side=1,at=c(1:nAges),labels=MeanLatAge,col="blue")

#  return(rep)
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
   par(mar=c(4,4,1,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   lm <- layout(matrix(c(1:3),nrow=3,byrow=TRUE))
   layout.show(lm)

   nice.ts.plot(MeanWatAge,MatAge,
        ylab="Natural Mortality (q)",xlab="Weight (kg)")
   legend("topleft",bty='n',cex=1.6,legend="M")

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
           ylab="Fishing Mortality (q)",xlab="Weight (kg)")
      points(c(3*kgperlb,10*kgperlb,15*kgperlb,20*kgperlb),c(0,0,0,0),col="red",pch='|',cex=2)
      legend("topleft",bty='n',cex=1.6,legend=paste("F(",r,")",sep=""))
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
   double.lines(csm.rep$mort.mat$wt,csm$est[qM],bcol="red",fcol="orange",lwd=7)
   sd.bars(csm.rep$mort.mat$wt,csm$est[qM],2.0*csm$std[qM],lwd=2,col="seagreen")

   maxF = 1.2*max(AveFatAge,csm$est[qF])
   print(maxF)
   nice.ts.plot(MeanWatAge,AveFatAge,bcol="black",fcol="lightgray",
                legend=c(" R2"," R4"), ylim = c(0,maxF),
                xlab="Weight (kg)", 
                ylab=substitute(paste("Fishing Mortality ",(q^{-1}))))
   double.lines(csm.rep$mort.mat$wt,csm$est[qF],bcol="red",fcol="orange",lwd=7)
   sd.bars(csm.rep$mort.mat$wt,csm$est[qF],2.0*csm$std[qF],col="seagreen")


   save.png.plot("csm_MF",width=width,height=height)
   par(old.par)
   
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
   }
}

##########
do.it.all=function()
{
   graphics.off()
   mfcl.ypr(region=2)->tmp
   mfcl.ypr(region=4)->tmp
   plot.mfcl.mortality()->tmp
}
