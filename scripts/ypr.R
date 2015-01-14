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

# from Dalzell
# The current minimum size for YFT for commercial sale in Hawaii is 3 lb.
# Is it possible that the Y/R analysis could look at a range of potential minimum sizes.
# I believe Dave Itano has suggested 10 lb, but maybe we would want to also check out 
# what would be the impact of 15 lb and 20 lb on the Y/R.

kgperlb = 0.453592

ypr<-function(epoch=NULL,region=2,astar=NULL)
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
#  print(Fmult)
   nf = length(Fmult)
   YPR = vector(length=nf)
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
      YPR[m] = sum(Y)
   }


   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(5,4,0,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   lm <- layout(matrix(c(1:2),nrow=2,byrow=TRUE))
   layout.show(lm)
   plot(MeanWatAge,AveFatAge[region,],las=1,type='b',ylim=c(0,max(AveFatAge[region,])),
         ylab="Fishing Mortality",xlab="Weight (kg)")
   points(c(3*kgperlb,10*kgperlb,15*kgperlb,20*kgperlb),c(0,0,0,0),col="red",pch='|')#,cex=2)
   if(astar < nAges)
      abline(v=MeanWatAge[astar],col="green",lty="dotdash")

   plot(Fmult,YPR,las=1,type='b',ylim=c(0,max(YPR)),
         xlab="F multiplier",ylab="Yield per Recruit (kg)")
   abline(v=1.0,col="red",lty="dotdash")

   par(old.par)
  return(as.data.frame(cbind(YPR,Fmult)))
}

plot.yield<-function(epoch=NULL,region=2)
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

ypr.step = function()
{
   N = prevN*exp(-Z)	
   B = prevB + prevN*W + N*prevW + 2.0*prevN*prevW
}

