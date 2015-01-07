
MHI.sim<-function(r = NULL, T12 = 0.01, T21 = 0.001,dt=1.0,fr=2,seta=c(0,0),
                  F.fcng=FALSE,T.fcng=TRUE, USE.LOGS = FALSE,
                  yr1=1952,yr2=2012,cfile="five_gears.dat") 
{
   ntime = (yr2 - yr1 + 1)*4
   F = compute.F(yr1=yr1,yr2=yr2,cfile=cfile)
   print(dim(F))
   ngear = nrow(F)
   sumF = colSums(F)
   print(paste(yr1,yr2,ntime,dim(F)[2],length(sumF)))
   F.summary = summary(sumF)
   F.msy = F.summary[5] # 3rd Qu.
   print(F.summary)

   if (is.null(r))
   {
      print(F.summary)
      r = 2.0*F.msy
      print(paste("r =",r))
   }  
   if (!F.fcng)
   {
      sumF = rep(F.msy,length(sumF))
   #  print(sumF)
   }

   ss = 1.0/dt 
   print(paste(ss,dt))

   region.biomass = as.matrix(read.table("total_biomass.dat"))
   immigrant.biomass = region.biomass[fr,]
   print(paste(fr,mean(immigrant.biomass)))
   # average region 2 biomass =  79953 
   # average region 4 biomass = 455982
   if (!T.fcng)
   {
      mib = mean(immigrant.biomass)
      immigrant.biomass = rep(mib,length(immigrant.biomass))
   #  print(immigrant.biomass)
   }
   K = 0.5*immigrant.biomass[1] 
   B.msy = K*(1.0-F.msy/K)

   eta<-matrix(0.0,nrow=2,ncol=ntime)
#  print(paste("seta_1 =",seta[1]))
   eta[1,] = rnorm(ntime,sd=seta[1])
#  print(eta[1,])
   eta[2,] = rnorm(ntime,sd=seta[2])
#  print(eta[2,])

   yield = matrix(0,nrow=ngear,ncol=ntime)
  
   print(paste("USE.LOGS =",USE.LOGS))
   if (USE.LOGS)
   {
   ##############################################
   x1=vector(length=ntime)
   x1[1] = log(0.9*K)
   x2=vector(length=ntime)
   x2[1] = log(0.1*K)
   for (i in 2:ntime)
   {
      nextx1 = x1[i-1]
      nextx2 = x2[i-1]
      for (k in 1:ss)
      {
         prevx1 = nextx1 #prevLogN1 = nextLogN1;
         prevx2 = nextx2 #prevLogN2 = nextLogN2;
         prevN1 = exp(prevx1) #prevN1 = exp(prevLogN1);
         prevN2 = exp(prevx2) #prevN2 = exp(prevLogN2);
      #  nextLogN1      += dt*(r*(1.0 - prevN1/K) - sumFg   - T12 - r*prevN2/K);
         nextx1 = nextx1 + dt*(r*(1.0 - prevN1/K) - sumF[i] - T12 - r*prevN2/K)

      #  nextLogN2      += dt*(r*(1.0 - prevN2/K) - sumFg  - T12 - r*prevN1/K + T21*immigrant_biomass(t)/prevN2);
         nextx2 = nextx2 + dt*(r*(1.0 - prevN2/K) - sumF[i] -T12 - r*prevN1/K + T21*immigrant.biomass[i]/prevN2)

      }
      x1[i] = nextx1
      x2[i] = nextx2

     # log_total_mean_pop = log( 0.5*(exp(pop11) + exp(pop21) + exp(pop21) + exp(pop22)) );

      total.mean.pop = 0.5*(exp(x1[i-1])+exp(x2[i-1])+exp(x1[i])+exp(x2[i]))
      for (g in 1:ngear)
      {
      #  log_pred_yield(g) =  ft(g) + log_total_mean_pop;
         yield[g,i] = F[g,i] * total.mean.pop
      }

   } # for (i in 2:ntime)


   N1 = exp(x1+eta[1,i])

   N2 = exp(x2+eta[2,i])
   }

   else
   {
   ##############################################
   N1=vector(length=ntime)
   N1[1] = 0.9*K
   N2=vector(length=ntime)
   N2[1] = 0.1*K
   for (g in 1:ngear)
      yield[g,1] = F[g,1]*(N1[1]+N2[1])
   
   for (i in 2:ntime)
   {
      prevN1 = N1[i-1]
      prevN2 = N2[i-1]
      for (k in 1:ss)
      {
         dN1 = dt*r*prevN1*(1.0 - prevN1/K) - dt*prevN1*(sumF[i]+T12) - r*dt*prevN1*prevN2/K
         dN2 = dt*r*prevN2*(1.0 - prevN2/K) - dt*prevN2*(sumF[i]+T12) - r*dt*prevN1*prevN2/K + dt*T21*immigrant.biomass[i]
         nextN1 = (prevN1+dN1)*exp(eta[1,i])
         nextN2 = (prevN2+dN2)*exp(eta[2,i])
         dtN = dt*(nextN1+nextN2)
         for (g in 1:ngear)
            yield[g,i] = yield[g,i] + F[g,i]*dtN
         
      }
      N1[i] = nextN1
      N2[i] = nextN2
   }
   ##############################################
   }

   sumY = colSums(yield)
   prop = vector(length=ntime)
   prop = N1/(N1+N2) #exp(x1)/(exp(x1)+exp(x2))
   print(summary(prop))
   tt <- seq((yr1+0.125),(yr2+0.875),0.25)
   print(length(tt))
   main = paste("T12=",T12,", T21=",T21," (",fr,")",
                ", S=(",seta[1],",",seta[2],")",sep="")
   print(main)
   sub = paste("r = ",r,", K= ",K,", dt = ",dt,", (ss = ",ss,")",sep="")
   print(sub)

   old.par <- par(no.readonly = TRUE) 
#  par(mar=c(5,4,4,7)+0.1)
#  par(mar=c(5,4,4,2)+0.1)


   lwd=5
   nice.ts.plot(tt,N1,bcol="blue",fcol="lightblue",lwd=lwd+2)
   double.lines(tt,N2,bcol="blue",fcol="lightblue",lwd=lwd+2,lty="dotted")
   double.lines(tt,sumY,bcol="darkred",fcol="red",lwd=lwd-2)
   title(main=main,sub=sub,ylab="Biomass (mt)",xlab="Year")
#  abline(h=B.msy,lwd=lwd,col="red")
   lines(range(tt),c(B.msy,B.msy),lwd=lwd,col="blue")
#  rug(tt)

   par("new"=TRUE)
   ytics = c(0.2,0.4,0.6,0.8,1.0)
   plot(range(tt),c(0,1),type='n',axes=FALSE,ann="FALSE")
   short.abline(c(range(tt)[1],par("usr")[2]),c(0.9,0.9),lwd=1,col="lightgreen")
   double.lines(tt,prop,bcol="darkgreen",fcol="lightgreen",lwd=lwd)
   axis(side=4,las=1,at=ytics,outer=FALSE,tcl=0.5,line=0,col.ticks="darkgreen")
   mtext("P1",side=4,line=0,col="darkgreen")
   abline(v=par("usr")[2],lwd=3,col="darkgreen")


   X11()
   lm = layout(matrix(c(1:ngear),ncol=1,byrow=TRUE))
   layout.show(lm)
   for (g in 1:ngear)
     nice.ts.plot(tt,yield[g,])

   dat.file = paste("simulated_catch.dat",sep="")
   print(paste("writing",dat.file))
   write(yield,file=dat.file,ncolumns=ntime)
   print(paste("finished",dat.file))

#  par(old.par)
#  return(as.matrix(cbind(tt,N1,N2,sumF,immigrant.biomass)))
   return(as.matrix(cbind(tt,t(yield))))
}

compute.F<-function(yr1=1952,yr2=2012,cfile="hdar_1952_2012.dat")
{
   eps.na = 1e-8
   obs.catch = read.table(file=cfile)
   nfish = nrow(obs.catch)
   ntime = ncol(obs.catch)

   #  get rid of the NAs and compute maxima
   max.catch = vector(length=nfish)
   for (i in 1:nfish)
   {
      w = which(is.na(obs.catch[i,]))
      obs.catch[i,w] = eps.na
      max.catch[i] = max(obs.catch[i,])
   }
   print(max.catch)
   max.catch = max.catch/sum(max.catch)
   print(max.catch)
   print(sum(max.catch))

   #  compute a multipation factor for the increase over previous observation   
   f.mult = matrix(1.0,nrow=nfish,ncol=ntime)
   for (i in 1:nfish)
   {
      for (j in 2:ntime)  
      {
         if (obs.catch[i,j-1] <= 0.0)
           f.mult[i,j] = 1.0
         else  
           f.mult[i,j] = obs.catch[i,j]/obs.catch[i,j-1]
         if (is.na(f.mult[i,j]))
            print(paste(f.mult[i,j],obs.catch[i,j-1],obs.catch[i,j]))
      }
   }

   # compute fishing mortality time series
   F = matrix(1.0,nrow=nfish,ncol=ntime)
   for (i in 1:nfish)
   {
      for (j in 2:ntime)  
      {
         F[i,j] = F[i,j-1]*f.mult[i,j]
      }
      # scaled tothe maximum catch
      F[i,] = max.catch[i]*F[i,]/max(F[i,])
   }

  
#  ret=list(obs.catch=obs.catch,f.mult=f.mult,F=F)
#  return(ret)
   return(F)
}

plot.region.biomass<-function(yr1=1952,yr2=2012,bfile="total_biomass.dat") 
{
   region.biomass = as.matrix(read.table(file=bfile))
   tt = seq((yr1+0.125),(yr2+0.875),0.25)
   print(length(tt))
   lwd=5

   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(3,4.5,0,0)+0.1)
   np = 2
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)

   nice.ts.plot(tt,region.biomass[2,],bcol="darkgreen",fcol="lightgreen",lwd=lwd,label="Region 2")
   dR2 = diff(region.biomass[2,])
   par("new"=TRUE)
   plot(tt[-1],dR2,type='l',lty="dotted",axes=FALSE,ann="FALSE")

   nice.ts.plot(tt,region.biomass[4,],bcol="darkgreen",fcol="lightgreen",lwd=lwd,label="Region 4")
   dR4 = diff(region.biomass[4,])
   par("new"=TRUE)
   plot(tt[-1],dR4,type='l',lty="dotted",axes=FALSE,ann="FALSE")

   print(summary(t(region.biomass[c(2,4),])))
   save.png.plot("MFCL_region_biomass",width=width,height=height)

   x11(width=width,height=height)
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   print(summary(t(rbind(dR2,dR4))))

   sd = sd(dR2)
   mean = mean(dR2)
   print(paste(mean,sd))
   breaks = 7 #seq(mean-3.25*sd,mean+3.25*sd,0.5*sd)
   print(breaks)
   hist(dR2,breaks=breaks,main="delta R2",freq=FALSE,las=1)
   x = seq(-3*sd,3*sd,0.1)
   lines(x,dnorm(x,mean=mean,sd=sd),col="blue")

   sd = sd(dR4)
   mean = mean(dR4)
   breaks = 7 #seq(-2.25*sd,2.25*sd,0.5*sd)
   hist(dR4,breaks=breaks,main="delta R2",freq=FALSE,las=1)
   x = seq(-3*sd,3*sd,0.1)
   lines(x,dnorm(x,mean=mean,sd=sd),col="blue")

   print(paste(sd(dR2),sd(dR4)))



   par(old.par)
#  return(region.biomass)
}

# compute initial estimates of sdlogF and sdlogYield as
# third quartile of first differneces in time series
std.comp<-function(yr1=1952,yr2=2012,cfile="hdar_1952_2012.dat")
{
   F = t(compute.F(yr1=yr1,yr2=yr2,cfile=cfile))
   wF = which(F <= 0)
   F[wF] = NA
   logF = log(F)
   dlogF = diff(logF)
   summF =summary(dlogF,na.rm=TRUE) 
#  print(summF[5,])

   Yield = t(read.table(file=cfile))
   wY = which(Yield <= 0)
   Yield[wY] = NA 
   logY = log(Yield)
   dlogY = diff(logY)
   summY = summary(dlogY,na.rm=TRUE)
#  print(summY[5,])
   return(list(sdlogF=summF[5,],sdlogYield=summY[5,]))
}

# logit(N1/(N1+N2)) = log(N1)-log(N2)
logit<-function(p)
{
   return(log(p/(1-p)))
}

alogit=function(alpha)
{
   return(1/(1+exp(-alpha)))
}

