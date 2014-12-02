require("R4MFCL")

ypr<-function(epoch=NULL,region=2)
{
   rep.file = paste(path.expand("~"),
                     "/Projects/HI-xfer/mfcl/yft/2014/basecase/plot-12.Rdata",
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

#  print(epoch)

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
