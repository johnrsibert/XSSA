require("mvtnorm")
gear.names = c("Tuna HL","Troll","Longline","Bottom/inshore HL","Aku Boat",
                "Klingon")
sgn = c("THL","Troll","LL","BHL","Aku")
have.xssams.R = TRUE
have.issams.R = FALSE
 
dLogN1=function(pN1,pN2,r, K, F, q, T12)
{
#  nextLogN1 += dt*(r*(1.0 - prevN1/K) - sumFg - T12 - 2.0*(1.0-q)*r*prevN2/K);
   dLN1 = r*(1.0 - pN1/K) - F - T12 - 2.0*(1.0-q)*r*pN2/K;
   return(dLN1)
}

dLogN2=function(pN1,pN2,r, K, F, q, T12, T21)
{
#  nextLogN2 += dt*(r*(1.0 - prevN2/K) - sumFg - T12 - 2.0*q*r*prevN1/K + T21*immigrant_biomass(t)/prevN2);
   dLN2 = r*(1.0 - pN2/K) - F - T12 - 2.0*q*r*pN1/K + T21/pN2;
   return(dLN2)
}



step=function(pN1,pN2,r, K, F, q, T12, T21, dt, s)
{
#  print(paste("p:",pN1,pN2))
   dLN1 = dLogN1(pN1,pN2,r, K, F, q, T12)
   dLN2 = dLogN2(pN1,pN2,r, K, F, q, T12, T21)
#  print(paste("d:",dLN1,dLN2))
   nN = vector(length=2,mode="numeric")
#  print(paste("s:",s[1],s[2]))
   nN[1] = exp(log(pN1) + dLN1 * dt + s[1])
   nN[2] = exp(log(pN2) + dLN2 * dt + s[2])
#  print(paste("nN:",nN[1],nN[2]))

#  nstep = 8
#  sdt = dt/nstep
#  print(paste("dt sdt:",dt,sdt))
#  tnN = vector(length=2,mode="numeric")
#  nextN1 = pN1
#  nextN2 = pN2
#  for (ss in 1:nstep)
#  {
#     prevN1 = nextN1
#     prevN2 = nextN2
#     dLN1 = dLogN1(prevN1,prevN2,r, K, F, q, T12)
#     dLN2 = dLogN2(prevN1,prevN2,r, K, F, q, T12, T21)
#     nextN1 = exp(log(prevN1) + dLN1 * sdt)
#     nextN2 = exp(log(prevN2) + dLN2 * sdt)
#  }
#  tnN[1] = nextN1*exp(s[1])
#  tnN[2] = nextN2*exp(s[2])
#  print(tnN)
#  print(paste("s:",s[1],s[2]))
#  print(tnN+s)
#  print(paste("tnN",tnN[1],tnN[2]))
   
   return(nN)
}

imp.step=function(pN1,pN2,r, K, F, q, T12, T21, dt, s)
{
   nN = vector(length=2,mode="numeric")
   nN[1] = pN1*(1.0+dt*r-dt*(F+T12)-dt*(1.0-q)*2.0*r*pN2/K)/
               (1.0+dt*r*pN1/K)
   nN[2] = (pN2*(1.0+dt*r-dt*(F+T12)-dt*q*2.0*r*pN1/K)+dt*T21)/
               (1.0+dt*r*pN2/K)
   nN = nN*exp(s)

#  print(paste("nN:",nN[1],nN[2]))
   return(nN)
}

daN1=function(pN1,pN2,r, K, F, q, T12)
{
   dN = pN1*(r*(1.0 - pN1/K) - F - T12) - 2.0*(1.0-q)*r*pN1*pN2/K;
   return(dN)
}

daN2=function(pN1,pN2,r, K, F, q, T12, T21)
{
   dN = pN2*(r*(1.0 - pN2/K) - F - T12) - 2.0*q*r*pN1*pN2/K + T21;
   return(dN)
}


stepa=function(pN1,pN2,r, K, F, q, T12, T21, dt, s)
{
#  print(paste("p:",pN1,pN2))
   dN1 = daN1(pN1,pN2,r, K, F, q, T12)
   dN2 = daN2(pN1,pN2,r, K, F, q, T12, T21)
#  print(paste("d:",dN1,dN2))

   nN = vector(length=2)
   nN[1] = pN1 + dN1 * dt
   nN[2] = pN2 + dN2 * dt
#  print(paste("n:",nN[1],nN[2]))
   return(nN)
}

zcompa=function(pN1,pN2,r, K, FF, q, T12, T21) 
{

   dN1 = daN1(pN1,pN2,r, K, FF, q, T12)
   dN2 = daN2(pN1,pN2,r, K, FF, q, T12, T21)
   z = dN2/dN1
#  print(paste(dN1,dN2,z))
   return(z)
}

# test reaterVis example
trast=function()
{
   proj <- CRS('+proj=longlat +datum=WGS84')
   df <- expand.grid(x=seq(-2, 2, .01), y=seq(-2, 2, .01))
   print(tail(df))
   df$z <- with(df, (3*x^2 + y)*exp(-x^2-y^2))
   print(tail(df))
   r1 <- rasterFromXYZ(df, crs=proj)
#  vectorplot(r1)
   vectorplot(r1,xlab="X",ylab="Y",colorkey=FALSE,countour=FALSE)
}

NNphase=function(r=0.3, K=1.0, FF = 0.007, q=0.54, T12=0.01, T21=0.002, dt = 1.0, s=c(0.0,0.0,0.0))
{
   USE_rasterVis = FALSE
   if (USE_rasterVis)
   {
   res = 0.01
   df =  expand.grid(pN1=seq(0.1, 1.1, res), pN2=seq(0.1, 1.1, res))
   print(tail(df))

#  df$z =  with(df, deparse(zcompa), r, K, FF, q, T12, T21)
   ngrid = nrow(df)
   print(names(df))
   z = vector(length=ngrid)
   for (i in 1:ngrid)
   {
      z[i] = zcompa(df[i,1],df[i,2],r, K, FF, q, T12, T21)
   }
   df = cbind(df,z)
   print(tail(df))

   proj <- CRS('+proj=longlat +datum=WGS84')
   rr = rasterFromXYZ(df,crs=proj)

   vectorplot(rr,xlab="N1",ylab="N2",main=paste("q =",q))
   }

   #### here
   else
   {
   require("plotrix")
   res = 0.05
   iN1 = seq(0.0, 1.1, res)
   iN2 = seq(0.0, 1.1, res) 
   nc = length(iN1)
   nr = length(iN2)
   prop = matrix(ncol = nc, nrow = nc)
   for (i in 1:nr)
   {
      for (j in 1:nc)
      {
         prop[i,j] = iN1[i]/(iN1[i] + iN2[j] + 1e-8)
      }
   }

   pN1 = seq(0.00, 1.1, res)
   pN2 = seq(0.00, 1.1, res) 
   df = expand.grid(pN1,pN2)
   ngrid = nrow(df)
   u = vector(length=ngrid)
   v = vector(length=ngrid)
   for (i in 1:ngrid)
   {
      u[i] = daN1(df[i,1],df[i,2],r, K, FF, q, T12) 
      v[i] = daN2(df[i,1],df[i,2],r, K, FF, q, T12, T21)
   }
#  plot(pN1,pN2,type='n',xlab="N1",ylab="N2")
   plot(c(0.0,1.1),c(0.0,1.15),type='n',
        xlab=expression("N"[1][,][1]),ylab=expression("N"[2][,][1]))
#  image(iN1,iN2,prop,zlim=c(0,1),col=heat.colors(5),add=TRUE)

#  pl = 0.9
#  p = pN1*(1.0/pl - 1.0)
   p = aPropL(pN1,0.9)
   double.lines(pN1,p,lwd=15,bcol="darkgreen",fcol="lightgreen")

#  contour(iN1,iN2,prop,zlim=c(0,1),col="darkgreen",add=TRUE)

   pl = seq(0,1,0.1)
   npl = length(pl)
   nn = length(pN1)
   for (i in 1:npl)
   {
      N2 = aPropL(pN1,pl[i])
      lines(pN1,N2,col="darkgreen",lwd=2,lty="dashed")
      text(pN1[nn],N2[nn],pl[i],col="darkgreen",adj=c(0.5,0.5))
   }
   lines(c(0.0,K),c(K,0.0), col="red",lwd=3, lty="dashed")
#  title(main=paste("r=",r,", K=", K, ", F=", FF, ", q=", q,
#                   ", T12=", T12, ",T21=", T21,sep=""),
#         line=-1,font.main=1,cex=0.8)
#  legend("topleft",bty='n',cex=1.6,legend=paste("T21/q=",T21/q,sep=""),text.font=2)
 
#  mat = cbind(u,v,df)
#  print(tail(mat))
   vectorField(u,v,xpos=df[,1],ypos=df[,2],vecspec="lonlat",
               headspan=0.05,scale=5,col="blue")
   }
}

aPropL = function(N1,p)
{
   N2 = N1*(1.0/(p+1e-8) - 1.0)
}

qcomp.ts=function(q)
{
   n = length(q)
#  width = 6.5
#  height = 9.0
#  x11(width=width,height=height)
#  old.par = par(no.readonly = TRUE) 
#  par(mar=c(4,4,0,0)+0.1)
#  lm = layout(matrix(c(1:n),ncol=1,nrow=length(q),byrow=TRUE))
   for(k in 1:length(q))
   {
      xssams.sim(q=q[k])
      title(main=paste("q = ",q[k],sep=""),
                          line=2,font.main=1,cex=0.8)
   }
}

qcomp.phase.4=function()
{
   q = c(0.25,0.75)
   T21 = c(0.02,0.0002)
   n = length(q)*length(T21)
   width = 9.0
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,4,0,0)+0.1)

   lm = layout(matrix(c(1:n),ncol=length(q),nrow=length(T21),byrow=TRUE))
   
   for (t in 1:length(T21))
   {
      for(k in 1:length(q))
      {
         NNphase(q=q[k],T21=T21[t])
         rat = round(T21[t]/q[k],3)
         title(main=paste("q=",q[k],", T21=",T21[t],", T21/q=",rat,sep=""),
                          line=-1,font.main=1,cex=0.8)
      }
   }
   save.png.plot(paste("qcomp_",n,sep=""),width=width,height=height)

   par(old.par)
}

qcomp.phase.2=function(T21=0.001)
{
   q = c(0.25,0.75)
   n = length(q)
   height = 6.5
   width = 2 * height
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,4,0,0)+0.1)

   lm = layout(matrix(c(1:n),ncol=length(q),nrow=1,byrow=FALSE))
   
      for(k in 1:length(q))
      {
         NNphase(q=q[k],T21=T21)
         rat = round(T21/q[k],5)
         title(main=paste("q=",q[k],", T21=",T21,", T21/q=",rat,sep=""),
                          line=-1,font.main=1,cex=0.8)
      }
   save.png.plot(paste("qcomp_",n,sep=""),width=width,height=height)

   par(old.par)
}

plot.NN.ts=function(pop,r=NULL,K,ib,save.graphics=TRUE)
{
   ntime = nrow(pop)

   lprop = pop[,1]/pop[,3]

   width = 9.0
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,5,0,4)+0.1)
   options(scipen=6)
   x = c(1:ntime)
   xrange=nice.ts.plot(x,pop[,1:3],legend=colnames(pop),xlab="t",ylab="N (mt)")
   abline(h=K,lwd=2,lty="dotdash",col="blue")
   lines(x,pop[,4],lwd=3,col="blue",lty="dashed")
   lines(x,pop[,5],lwd=3,col="blue",lty="dashed")
   text(x[ntime],K," K",adj=c(0,0),col="blue")
   #if (!is.null(r))
   #   title(main=paste("r = ",r,sep=""),line=-1)

   par("new"=TRUE)
   plot(x,lprop,lwd=3,type='l',col="red",ylim=c(0,1),
        ann=FALSE,axes=FALSE,xlim=xrange)
   text(x[ntime],lprop[ntime]," p",adj=c(0,0),col="red")
   abline(h=0.9,lwd=2,lty="dotdash",col="red")
   axis(4,col="red",ylab="p",col.axis="red")
   #abline(v=par("usr")[2],lwd=2,col="red")
   mtext("p",side=4,col="red",line=0.1)

   par("new"=TRUE)
#  plot(x,ib,type='n',ann=FALSE,axes=FALSE,
#       ylim=c(0.0,max(ib)),xlim=xrange)
#  double.lines(x,ib,bcol="purple4",fcol="purple1",lwd=5)
   plot(x,ib,lwd=3,type='l',col="purple1", 
        ylim=c(0.0,max(ib)), ann=FALSE,axes=FALSE,xlim=xrange)
   text(x[ntime],ib[ntime]," T21",adj=c(0,0),col="purple1")
   axis(4,line=-2,col="purple1",col.axis="purple1")

   if (save.graphics)
      save.png.plot("NNts",width=width,height=height)
   par(old.par)
}

plotN1N2=function(ntime=100,r=0.3, K=1.0, F = 0.007, q=0.54, 
                  T12=0.01, T21=0.002, p=0.9, 
                  dt = 1.0, s=c(0.0,0.0,0.0))
{
   pop = matrix(ncol=3,nrow=ntime)
   colnames(pop)=c("N1","N2","N1+N2")
#  colnames(pop)=c(expression("N"[1]),
#                  expression("N"[2]),
#                  expression("N"[1]+"N"[2]))
   pop[1,1] = p*K
   pop[1,2] = (1-p)*K
   pop[1,3] = pop[1,1]+pop[1,2]
#  print(head(pop))
   for (t in 2:ntime)
   {
      N1 = pop[t-1,1]
      N2 = pop[t-1,2]
      tN = step(N1,N2,r,K,F,q,T12,T21,dt,s)
      pop[t,1:2] = tN
      pop[t,3] = pop[t,1]+pop[t,2]
   #  print(pop[t,])
   }
  
   plot.NN.ts(pop,r,K)

}

plotN1N2a=function(ntime=100,r=0.3, K=1.0, F = 0.007, q=0.54, 
                  T12=0.01, T21=0.002, p=0.9, 
                  dt = 1.0, s=c(0.0,0.0,0.0))
{
   pop = matrix(ncol=3,nrow=ntime)
   colnames(pop)=c("N1","N2","N1+N2")
   pop[1,1] = p*K
   pop[1,2] = (1-p)*K
   pop[1,3] = pop[1,1]+pop[2,1]
#  print(head(pop))
   for (t in 2:ntime)
   {
      N1 = pop[t-1,1]
      N2 = pop[t-1,2]
      tN = stepa(N1,N2,r,K,F,q,T12,T21,dt,s)
      pop[t,1:2] = tN
      pop[t,3] = pop[t,1]+pop[t,2]
   #  print(pop[t,])
   }
   plot.NN.ts(pop,r,K)
##   title("Arithmetic")
}

a.log.disp=function(ntime=100,r=0.3, K=1.0, F = 0.007, q=0.54, 
                  T12=0.01, T21=0.002, p=0.9, 
                  dt = 1.0, s=c(0.0,0.0,0.0))
{

   plotN1N2(ntime, r, K, F, q, T12, T21, p, dt, s)
   x11()
   plotN1N2a(ntime, r, K, F, q, T12, T21, p, dt, s)

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


plot.propL.prior=function(propL.prior=0.9)
{
   p = seq(0.01,0.99,0.01)
   p = c(0.0001,p,0.9999)
   wp = which(p == propL.prior)

   Lp = logit(p)
   LpropL.prior = logit(propL.prior)

#  lm = layout(matrix(c(1:2),ncol=1,byrow=TRUE))
#  layout.show(lm)

   
   
   plot(c(0.0,1.0),c(0.0,2.0),type='n',xlab="x",ylab="p(L(x))")
   abline(v=propL.prior,col="red",lty="dotdash")
   Lsd = c(0.501,0.6,0.7,0.8,0.9,0.99,0.999)
   print(Lsd)
   asd = logit(Lsd)
   for (sd in asd)
   {
      pp = dnorm(Lp,mean=LpropL.prior,sd=sd)
   #  plot(Lp,pp)
      print(paste("sd",sd,alogit(sd)))
      double.lines(alogit(Lp),pp,lwd=5,fcol="lightblue",bcol="blue")
      text(propL.prior,pp[wp],alogit(sd),col="blue")
   }

   save.png.plot("propL_prior",width=6.5,height=6.5)
}




#dnorm(x, mean = 0, sd = 1, log = FALSE)
# standard Cauchy density
dtdist=function(x, mean = 0, sd = 1)
{
  xx = (x-mean)/sd
# print(paste(xx,xx*xx,(1.0+xx*xx)))
  dt = 1.0/(pi*(1.0+xx*xx))
  return(dt)
}


#  dvariable peN2 = 0.5*(log(TWO_M_PI*varlogPop(2)) + square(p22-nextLogN2)/varlogPop(2));
mydnorm=function(x,mean=0,sd=1)
{
   y = 0.5*(log(2*pi*sd*sd) + ((x-mean)/sd)^2)
   return(exp(-y))
}


plot.t.norm.mix=function(mean = 0,sd=1)
{
   x = seq(mean-4*sd,mean+4*sd,sd/10)
#  print(x)
   yn = dnorm(x,mean=mean,sd=sd)
   myn = mydnorm(x,mean=mean,sd=sd)
   print(paste("sum normal:",sum(yn)))
   intn = integrate(function(x) dnorm(x,mean=mean,sd=sd),lower=-Inf,upper=Inf)
   print(paste("  integral:",intn$value))
   print(paste("sum mynormal:",sum(myn)))
   intmyn = integrate(function(x) mydnorm(x,mean=mean,sd=sd),lower=-Inf,upper=Inf)
   print(paste("  integral:",intmyn$value))
   plot(x,myn,type='l')
   points(x,yn)
   tn = dtdist(x,mean=mean,sd=sd)
   print(paste("     sum t:",sum(tn)))
   intt = integrate(function(x) dtdist(x,mean=mean,sd=sd),lower=-Inf,upper=Inf)
   print(paste("  integral:",intt$value))
   lines(x,tn,col="red")

   pf = seq(0.1,0.9,0.1)
   for (pfat in pf)
   {
     robust=function() {pfat*tn+(1-pfat)*yn}
     lines(x,robust(),col="blue")
   }
#  lines(x,(pfat*tn+(1-pfat)*yn),col="blue")

}



# from MHI_sim.R
#compute.F<-function(yr1=1952,yr2=2012,cfile="../HDAR/hdar_1952_2012.dat",plot=TRUE)
compute.F<-function(obs.catch,plot=TRUE)
{
   eps.na = 1e-8
#  obs.catch = read.table(file=cfile)
   ngear = nrow(obs.catch)
   ntime = ncol(obs.catch)

   #  get rid of the NAs and compute maxima
   max.catch = vector(length=ngear)
   for (i in 1:ngear)
   {
      w = which(is.na(obs.catch[i,]))
      obs.catch[i,w] = eps.na
      max.catch[i] = max(obs.catch[i,])
   }
   max.catch = max.catch/sum(max.catch)

   #  compute a multipation factor for the increase over previous observation   
   f.mult = matrix(1.0,nrow=ngear,ncol=ntime)
   for (i in 1:ngear)
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
   F = matrix(1.0,nrow=ngear,ncol=ntime)
   for (i in 1:ngear)
   {
      for (j in 2:ntime)  
      {
         F[i,j] = F[i,j-1]*f.mult[i,j]
      }
      # scaled to the maximum catch
      F[i,] = max.catch[i]*F[i,]/max(F[i,])
   }

   if (plot)
   {
      x = c(1:ntime)
      width = 9.0
      height = 11.0
      x11(width=width,height=height)
      old.par = par(no.readonly = TRUE) 
      par(mar=c(3.5,4.5,0,0)+0.1)
      np = ngear
      lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
      layout.show(lm)
      for (g in 1:ngear)
      {
         if (g < ngear)
            nice.ts.plot(x,F[g,],bcol="orange4",fcol="orange",lwd=3,ylab="F")
         else
            nice.ts.plot(x,F[g,],bcol="orange4",fcol="orange",lwd=3,ylab="F",
                     xlab="Time")
         title(main=gear.names[g],line=-1.5)
      } 
      par(old.par)
   }
#  ret=list(obs.catch=obs.catch,f.mult=f.mult,F=F)
#  return(ret)
   return(F)
}

xssams.sim=function(r=1.2, K=200000, q=0.54, T12=0.01, T21=0.002,
                  dt=1.0, sN=c(0.0,0.0,0.0), fr=2, p=0.9, 
                  F=NULL, sC=c(0.0,0.0,0.0,0.0,0.0), F.mult=1.0,
                  do.plot=TRUE, save.graphics=FALSE, do.est=FALSE,
                  seed=NULL)
{
   if (!is.null(seed))
   {
      set.seed(seed)
      print(paste("random number seed set to",seed))
   }

   print(paste("r = ",r,", dt =",dt,sep=""))
   if (dt == 0.25)
   {
      obs.catch = t(read.table(file="../run/five_gears_q_1952_2012.dat"))
      biomass.file ="../run/total_biomass_q.dat"
      region.biomass = as.matrix(read.table(biomass.file))
   }
   else if (dt == 1.0)
   {
      obs.catch = t(read.table(file="../run/five_gears_a_1952_2012.dat"))
      biomass.file = "../run/total_biomass_a.dat"
      region.biomass = as.matrix(read.table(biomass.file))
   }
   else
   {
      print(paste("dt = ",dt," Unknown.",sep=""))
      return(FALSE)
   }
#  print("region.biomass:")
#  print(dim(region.biomass))
#  print(head(t(region.biomass)))
   immigrant.biomass = T21*region.biomass[fr,]

   if (is.null(F))
   {
      print("computing F")
      F.scale=0.015 # observed and predicted catch similar
      F.matrix = F.scale*compute.F(t(obs.catch),plot=FALSE) 
   }
   else
      F.matrix = F

   ntime = ncol(F.matrix)
   ngear = nrow(F.matrix)
#  print(paste(ngear,ntime))
   print(paste("F.mult =",F.mult))
   fm = seq(1, F.mult, (F.mult-1)/(ntime-1))
#  print(length(fm))
#  print(fm)
   
#  mean.F = mean(F.matrix,na.rm=TRUE)
#  print(paste("mean F = ",mean.F,",log(mean F) = ",log(mean.F),sep=""))
#  log.F = log(F.matrix+1)
#  print(paste("mean log (F) =",mean(log.F)))
   for (g in 1:ngear)
   {
      F.matrix[g,] = F.matrix[g,]*fm
   }
   mean.F = mean(F.matrix,na.rm=TRUE)
   print(paste("mean F = ",mean.F,",log(mean F) = ",log(mean.F),sep=""))
   log.F = log(F.matrix+1)
   print(paste("mean log (F) =",mean(log.F)))

   sum.F = colSums(F.matrix)
#  print(sN)
   cov = matrix(nrow=2,ncol=2)
   cov[1,1] = sN[1]
   cov[2,2] = sN[2]
   cov[1,2] = sN[3]*cov[1,1]*cov[2,2]
   cov[2,1] = cov[1,2]
#  print(cov)
   s <- rmvnorm(n=ntime, sigma=cov)
#  print(dim(s))
#  print(colMeans(s))
#  print(var(s))

   pop = matrix(ncol=5,nrow=ntime)
   colnames(pop)=c(" N1"," N2"," N1+N2","tN1","tN2")
   pop[1,1] = p*K*exp(s[1,1])
   pop[1,2] = (1-p)*K*exp(s[1,2])
   pop[1,3] = pop[1,1]+pop[1,2]
   pop[1,4] = pop[1,1]
   pop[1,5] = pop[1,2]
   for (t in 2:ntime)
   {
   #  print(paste("t =",t))
      N1 = pop[t-1,1]
      N2 = pop[t-1,2]
      
      tN = step(N1,N2,r,K,sum.F[t],q,T12,immigrant.biomass[t],dt,s[t,])
      pop[t,1:2] = tN[1:2]
      pop[t,3] = pop[t,1]+pop[t,2]

      N1 = pop[t-1,4]
      N2 = pop[t-1,5]
      tN = imp.step(N1,N2,r,K,sum.F[t],q,T12,immigrant.biomass[t],dt,s[t,])
      
      pop[t,4:5] = tN[1:2]
   }
  
   sC.mat = matrix(ncol=ngear,nrow=ntime)
   for (g in 1:ngear)
   {
      sC.mat[,g] = rnorm(ntime,mean=0.0,sd = sC[g])
   }
#  print("sC.mat:")
#  print(dim(sC.mat))
#  print(head(sC.mat))
#  print(sC.mat[1,])

   pred.catch = matrix(ncol=ngear,nrow=ntime)
   for (t in 1:ntime)
   {
      pred.catch[t,] = round(obs(pop[t,3],F.matrix[,t],sC.mat[t,]))
   }

#  print("pop:")
#  print(dim(pop))
#  print(head(pop))
#  print("obs.catch:")
#  print(dim(obs.catch))
#  print(head(obs.catch))
#  print("pred.catch:")
#  print(dim(pred.catch))
#  print(head(pred.catch))

   if (do.plot)
   {
      plot.NN.ts(pop,r,K,immigrant.biomass,save.graphics)
      title=paste("r=",r,", K=", K, ", Fx", F.mult, ", q=", q,
                  ", T12=", T12, ",T21=", T21,",sN=(", sN[1],",",sN[2],")",
                  ", seed=",seed,
                  sep="")
      title(main=title,cex.main=0.9,line=2.5)#sub=title)
      plot.catch.ts(obs.catch,pred.catch,save.graphics)
   }

   if (do.est)
   {
      dfile = "../run/xssams.dat"
      print(dfile)
   
      cat.number=function(v)
      {
         cat(paste(" ",v,"\n",sep=""),file=dfile,append=TRUE)
      }
      cat.string=function(s)
      {
         cat(paste(s,"\n",sep=""),file=dfile,append=TRUE)
      }
      cat.vector=function(v)
      {
         for (i in 1:length(v))
         {
            cat(paste(" ",v[i],sep=""),file=dfile,append=TRUE)
         }
         cat("\n",file=dfile,append=TRUE)
      }
      cat.matrix=function(m)
      {
          for (i in 1:nrow(m))
             cat.vector(m[i,])
      }
      cat("# xssams simulation output:\n",file=dfile)
      cat.string("#")
      cat.string("# Number of fishing gears")
      cat.number(ngear)
      cat.string("# Number of time periods")
      cat.number(ntime)
      cat.string("# Time step (years)")
      cat.number(dt)
      cat.string("# Catch data")
      cat.string(paste("# simulated catch; F.mult = ",F.mult,sep=""))
      cat.matrix(t(pred.catch))
      cat.string(paste("#",biomass.file))
      cat.matrix(region.biomass)
      cat.string("# MFCL forcing region number")
      cat.number(2)
      cat.string("# use mean forcing")
      cat.number(0)
      cat.string("#")
      cat.string("# T12  phase  initial value")
      cat.vector(c(2,0.01))
      cat.string("# T21  phase  initial value")
      cat.vector(c(2,0.02))
      cat.string("# r    phase  initial value")
      cat.vector(c(2,r))
      cat.string("# K    phase  initial value")
      cat.vector(c(2,K))
      cat.string("# sdlogF    phase  initial values")
      cat.vector(c(1,0.3)) #, 0.3, 0.3, 0.3, 0.3))
      cat.string("# sdlogPop    phase  initial values")
      cat.vector(c(2,0.1)) #,sN[2]+0.001))
   #  cat.string("# rho phase initial value")
   #  cat.vector(c(-1,sN[3]))
      cat.string("# sdlogYield  phase initial values")
      cat.vector(c(1,1.0))
      cat.string("# meanProportion_Local phase initial value")
      cat.vector(c(-1,p))
      cat.string("# sdLproportion_local phase initial value")
      cat.vector(c(-1,    1.3863))
      cat.string("# qProp phase initial value")
      cat.vector(c(2,q))
      cat.string("# robust yield likelihood")
      cat.string("# use pfat_phase pfat initial values")
      cat.vector(c(1,   -1,         0.07, 0.07, 0.07, 0.07, 0.07))
  } 

}

obs = function(NN,F,s)
{
   ngear = length(F)
   pred.catch = vector(length=ngear)
   for (g in 1:ngear)
   {
      lpc = log(NN) + log(F[g])  + s[g]
      pred.catch[g] = exp(lpc)
   }
   return(pred.catch)
}


plot.catch.ts=function(obs.catch,pred.catch,save.graphics=TRUE)
{
   ntime = nrow(pred.catch)
   ngear = ncol(pred.catch)
   x = c(1:ntime)
   width = 9.0
   height = 11.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,5,0,1)+0.1)
   np = ngear+1
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   for (g in 1:np)
   {
      if (g < np)
      {
         nice.ts.plot(x,pred.catch[,g],xlab="",ylab="Catch (mt)",
                  bcol="darkgreen",fcol="lightgreen",lwd=3)
         points(x,obs.catch[,g],col= "darkgreen",pch=16)
         title(main=gear.names[g],line=-1)
      }
      else
      {
         y = rowSums(pred.catch)
         nice.ts.plot(x,y,xlab="t",ylab="",
                  bcol="darkgreen",fcol="lightgreen",lwd=3)
         y = rowSums(obs.catch,na.rm=TRUE)
         points(x,y,col= "darkgreen",pch=16)
         title(main="Total",line=-1)
      }
   } 

   if (save.graphics)
      save.png.plot("catchts",width=width,height=height)
}

read.rep=function(file,ntime,ngear,dt)
{
   print(paste("Scanning file",file))
   rep = scan(file,what="character")

   p = grep("logT21:",rep)
   logT12 = as.numeric(rep[p+1])
   p = grep("logT12:",rep)
   logT21 = as.numeric(rep[p+1])
   p = grep("logr:" ,rep)                 
   logr = as.numeric(rep[p+1])
   p = grep("logK:",rep)
   logK = as.numeric(rep[p+1])
   p = grep("LmeanProportion_local:" ,rep)
   Lpropl = as.numeric(rep[p+1])
   p = grep("logsdLProportion_local:",rep)
   logsdLpropl = as.numeric(rep[p+1])
   p = grep("prop:" ,rep)                 
   prop = as.numeric(rep[p+1])
   p = grep("logsdlogF:",rep)
   logsdlogF = as.numeric(rep[p[1]+1])
   p = grep("logsdlogPop:" ,rep)          
   logsdlogPop = as.numeric(c(rep[p[1]+1],rep[p[1]+2]))
   p = grep("logsdlogYield:" ,rep)        
   logsdlogYield = as.numeric(rep[p[1]+1])
   p = grep("qProp:",rep)
   qProp = as.numeric(rep[p+1])

   res = grep("Residuals:",rep)
   print(res)
   afters = grep("after",rep)
   print(afters)

   max.counter = length(res)
   counter = max.counter
   print(paste(max.counter, "blocks found:"))

   ncol = (3*ngear+6)
   diag = matrix(nrow=ntime,ncol=ncol)
   cnames = vector(length=ncol)

   fc = res
   for (c in 1:ncol)
   {
      fc = fc + 1
      #  print(paste(fc,log[fc])) 
      cnames[c] = rep[fc]
   }
   print(cnames)
   for (t in 1:ntime)
   {
      for (c in 1:ncol)
      {
         fc = fc + 1
         diag[t,c] = as.numeric(rep[fc],ngear)
      }
   }
   colnames(diag)=cnames

   return(list(logT12=logT12,
      logT21=logT21,
      logr=logr,
      logK=logK,
      Lpropl=Lpropl,
      logsdLpropl=logsdLpropl,
      prop=prop,
      logsdlogF=logsdlogF,
      logsdlogPop=logsdlogPop,
      logsdlogYield=logsdlogYield,
      qProp=qProp,
      diag=as.data.frame(diag)))

}


xssams.rep.sim=function(rep.file,ntime,ngear,dt)
#xssams.sim=function(r=0.12, K=200000, q=0.54, T12=0.001, T21=0.0002,
#                  dt=1.0, sN=c(0.0,0.0,0.0), fr=2, p=0.9, 
#                  F=NULL, sC=c(0.1,0.1,0.1,0.1,0.1), F.mult=0.0025,
#                  do.plot=TRUE, save.graphics=FALSE, do.est=FALSE)
{
   rep=read.rep(rep.file,ntime,ngear,dt)
   xssams.sim(r=exp(rep$logr), K=exp(rep$logK), q=rep$qProp,
               T12=exp(rep$logT12), T21=exp(rep$logT21),
               dt=dt, sN=c(exp(rep$logsdlogPop),0.0), fr=2, p=0.9, 
               F=exp(rep$diag[,7:11]), sC=c(0.1,0.1,0.1,0.1,0.1), 
               F.mult=1.0,do.plot=TRUE, save.graphics=FALSE, do.est=FALSE)
}


b27.test.sim=function(
# Status block 27
# nll = -404.11
# nvar = 427
# current phase = 2
#  logT12 = -4.6211 (1)
   T12 = 0.0098419,
#  logT21 = -4.8741 (1)
   T21 = 0.007642,
# logr = 0.84139 (1)
#    r = 2.3196,
     r=2.1,
# logK = 10.343 (1)
     K = 31039,
#     logsdlogF: -0.073259 (1)
#        sdlogF: 0.92936
#   logsdlogPop:  -3.3826 -3.6009 (1)
#      sdlogPop:  0.033958 0.027299
# rho = 0 (0)
# logsdlogYield: -2.7842 (1)
#    sdlogYield: 0.061778
# LmeanProportion_local = 2.1972 (0)
#                  prop = 0.9
# logsdLProportion_local = 1.0006 (0)
#    sdLProportion_local = 2.72
# pfat =  0.07 0.07 0.07 0.07 0.07 (0)
 q = 0.59224,
                  dt=1.0, 
                  sN=c(0.033958,0.027299,0),
                  fr=2, 
                  p=0.9, 
                  F=NULL, 
                  sC=rep(0.061778,5),
                  F.mult=0.015,
                  do.plot=TRUE, 
                  save.graphics=FALSE, 
                  do.est=FALSE)
{
   xssams.sim(r, K, q, T12, T21, dt, sN, fr, p, F, sC, 
               F.mult, do.plot, save.graphics, do.est)

}

#plot.error=function(x,y,sd,bcol,fcol,mult=2)
#{
#   if (capabilities("cairo"))
#   {
#      sdyu = exp(log(y)+mult*sd)
#      sdyl = exp(log(y)-mult*sd)
#      frgb = col2rgb(fcol)/255
#      polygon(c(x,rev(x)),c(sdyl,rev(sdyu)),
#              border=bcol,lty="dashed",lwd=1,
#              col=rgb(frgb[1],frgb[2],frgb[3],0.5))
#   }
#   else
#      polygon(c(x,rev(x)),c(sdyl,rev(sdyu)),
#              border=bcol,lty="dashed",lwd=1,col=fcol)
#
#}


plot.diagnostics=function(dat=NULL,file="diagnostics.dat",dt,ngear,
                 sdlogPop, sdlogYield, sdlogF, sdlogQ, K, r, T12, q,
                 plot.Fmort,plot.prod, plot.impact, plot.Discr,
                 devices,block)
{
   print(paste("Block: ",block))
   print(is.data.frame(dat))
   print(head(dat))
   if (is.null(dat))
     dat = read.table(file=file,header=TRUE)

   start.year = 1952
   print(start.year)
   if (dat$t[1] == 1)
   {
      dat$t = (start.year-0.4*dt +  dat$t*dt)
   #  print(names(dat))
   }

   wy5 = which((floor(dat$t%%5))==0)
   ncol = ncol(dat)
#  print("Names of all variables:")
#  print(names(dat))
   gear.col = 6
   dd = c(2:3,(gear.col+1):ncol)
#  print("Names of log transformed variables:")
#  print(names(dat)[dd])
   dat[,dd] = exp(dat[,dd])

   #gear.names = c("TunaHL","Troll","Longline","Bottom/inshore HL","AkuBoat")
   title.line = -1
   lwd = 3
   #sd.lwd = 3
   #sd.lty = "dotted"
   #old.par = par(no.readonly = TRUE) 

   F.ndx=grep("F",names(dat)) 
   predC.ndx=grep("predC",names(dat)) 
   obsC.ndx=grep("obsC",names(dat)) 

   # biomass plots
   d = 1
   xpos = 0
   ypos = 0
   if (devices[d] > 0)
   {
      s = dev.set(devices[d])
   }
   else
   {
      width = 9.0
      height = 9.0
      xpos = 100
      ypos = 0
      x11(width=width,height=height,xpos=xpos,ypos=ypos,title="Biomass")
      devices[d] = dev.cur()
   }

   ntime = length(dat$t)
   y = matrix(0.0,nrow=ntime,ncol=3)
   y[,1] = dat$pop1
   y[,2] = dat$pop2
   y[,3] = dat$pop1 + dat$pop2
   plot.biomass(dat$t,y,sd=sdlogPop,block=block,K=dat$K,
                propL=dat$propL,forcing=dat$forcing)


#   par(mar=c(4,5,0,4)+0.1)
#   ntime = length(dat$t)
#   x = dat$t
#   y = matrix(nrow=ntime,ncol=3)
#   y[,1] = dat$pop1
#   y[,2] = dat$pop2
#   y[,3] = dat$pop1 + dat$pop2
# # legend = c(" N1"," N2"," N1+N2")
#   legend = c(parse(text=paste("N","[1]")),
#              parse(text=paste("N","[2]")),
#              parse(text=paste("N","[1]","+N","[2]")))
#
#   options(scipen=6)
#   xrange=nice.ts.plot(x,y,ylab="Biomass (mt)")
#
#   lines(x,dat$K,lwd=2,lty="dotdash",col="blue",xlim=xrange)
#   text(x[ntime],dat$K[ntime]," K",adj=c(0,0.5),col="blue")
#
#   sdlogNN = sqrt(4.0*sdlogPop*sdlogPop) # is is probably not correct
#   plot.error(dat$t,y[,3],sdlogNN, 
#                    bcol="blue",fcol="lightblue")
#   plot.error(dat$t,y[,1],sdlogPop,
#                    bcol="blue",fcol="lightblue")
#   plot.error(dat$t,y[,2],sdlogPop,
#                    bcol="blue",fcol="lightblue")
#   lines(dat$t,y[,1],col="blue",lwd=5)
#   text(x[ntime],y[ntime,1],parse(text=paste("N","[1]")),adj=c(0,0),col="blue")
#   lines(dat$t,y[,2],col="blue",lwd=5)
#   text(x[ntime],y[ntime,2],parse(text=paste("N","[2]")),adj=c(0,0),col="blue")
#   lines(dat$t,y[,3],col="blue",lwd=5)
#   par("new"=TRUE)
#   plot(x,dat$propL,lwd=3,type='l',col="red",ylim=c(0,1),
#        ann=FALSE,axes=FALSE,xlim=xrange)
#   text(x[ntime],dat$propL[ntime]," p",adj=c(0,0),col="red")
#   abline(h=0.9,lwd=2,lty="dotdash",col="red")
#   axis(4,col="red",ylab="p",col.axis="red")
#   mtext("p",side=4,col="red",line=0.1)
#
#   par("new"=TRUE)
#   tT21 = parse(text=paste("T","[21]"))
#   plot(x,dat$forcing,lwd=3,type='l',col="purple", 
#        ylim=c(0.0,max(dat$forcing)), ann=FALSE,axes=FALSE,xlim=xrange)
#   text(x[ntime],dat$forcing[ntime],tT21,adj=c(0,0),col="purple")
#   axis(4,line=-2,col="purple",col.axis="purple")
##  axis(2,col="purple",col.axis="purple",pos=xrange[2])
#   mtext(tT21, side=4,col="purple",line=-1.5)
#   show.block.number(block,dat$t[1])

   # catch plots
   d = 2
   if (devices[d] > 0)
   {
      s = dev.set(devices[d])
   }
   else
   {
     width = 9.0
     height =11.0
     xpos = xpos + 50
     ypos = ypos + 50
     x11(width=width,height=height,xpos=xpos,ypos=ypos,title="Catch")
     devices[d] = dev.cur()
   }

   plot.catches(dat$t,dat[,c(gear.col+2*ngear+1):(gear.col+3*ngear)],
                      dat[,c(gear.col+ngear+1)  :(gear.col+2*ngear)],
                      sdlogYield,block) 
                      
#   par(mar=c(3,4.5,0,0)+0.1)
#   np = ngear+1
#   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
#   layout.show(lm)
#   sum.obs = vector(length=ntime,mode="numeric")
#   sum.pred = vector(length=ntime,mode="numeric")
#   for (g in 1:ngear)
#   {
#      pred.ndx = gear.col+  ngear+g
#      obs.ndx  = gear.col+2*ngear+g
#
#      nice.ts.plot(dat$t,dat[,pred.ndx],
#         bcol="darkgreen",fcol="lightgreen",lwd=lwd,ylab="Catch (mt)")
#      plot.error(dat$t,dat[,pred.ndx],sdlogYield,
#                 bcol="darkgreen",fcol="lightgreen")
#      lines(dat$t,dat[,pred.ndx],col="darkgreen",lwd=lwd+2)
#      points(dat$t,dat[,obs.ndx],col= "darkgreen",pch=3,cex=3)
#      title(main=gear.names[g],line=title.line)
#
#      sum.obs  = sum.obs  + dat[,obs.ndx]
#      sum.pred = sum.pred + dat[,pred.ndx]
#   }
#
#   nice.ts.plot(dat$t,sum.pred,
#                 bcol="darkgreen",fcol="lightgreen",lwd=lwd,ylab="Catch (mt)")
#   lines(dat$t,sum.pred,col="darkgreen",lwd=lwd+2)
#   points(dat$t,sum.obs,col= "darkgreen",pch=3,cex=3)
#   title(main="All Fleets",line=title.line)
#
#   show.block.number(block,dat$t[1],line=2)

#  plot.Fmort = TRUE
   if (plot.Fmort)
   {
       d = d + 1
       if (devices[d] > 0)
       {
          s = dev.set(devices[d])
       }
       else
       {
          width = 9.0
          height =11.0
          xpos = xpos + 50
          ypos = ypos + 50
          x11(width=width,height=height,xpos=xpos,ypos=ypos,title="Fishing Mortality")
          devices[d] = dev.cur()
       }
       par(mar=c(3,4.5,0,0)+0.1)
       np = ngear
       lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
       layout.show(lm)
       for (g in 1:ngear)
       {
       #  labelsY=parse(text=paste(seq(50,100,10), "^o ", "*N", sep="")) 
          Flab = parse(text=paste("F~(","y^-1",")",sep=""))
          nice.ts.plot(dat$t,dat[,(gear.col+g)], ylab=Flab,
                bcol="orange4",fcol="orange", lwd=lwd)

          plot.error(dat$t,dat[,(gear.col+g)],sdlogF,
                    bcol="orange4",fcol="orange")
          lines(dat$t,dat[,(gear.col+g)],col="orange4",lwd=lwd+2)
          title(main=gear.names[g],line=title.line)
       }
       show.block.number(block,dat$t[1],line=2)
   } #if (plot.Fmort)

#  production plot
   if (plot.prod)
   {
       d = d + 1
       if (devices[d] > 0)
       {
          s = dev.set(devices[d])
       }
       else
       {
          width = 9.0
          height =9.0
          xpos = xpos + 50
          ypos = ypos + 50
          x11(width=width,height=height,xpos=xpos,ypos=ypos,
               title="Production")
          devices[d] = dev.cur()
       }


       Fmort = rowSums(dat[,F.ndx])
       obsC = rowSums(dat[,obsC.ndx])
       predC = rowSums(dat[,predC.ndx])
       plot.production(Fmort,obsC,predC,dat$t, r,K,block)

#       Fmort = rowSums(dat[,F.ndx])
#       F.max=max(Fmort,na.rm=TRUE)
#       if (F.max > r)
#         F.max = max(F.max,r)
#       Fyield = seq(0,F.max,0.01*F.max)
#       yield = Fyield*K*(1.0-Fyield/r) # equilibirum yield at F
#       obsC = rowSums(dat[,obsC.ndx])
#       predC = rowSums(dat[,predC.ndx])
#       xrange = c(0,F.max)
#       yrange = c(0,max(obsC,predC,yield,na.rm=TRUE))
#       Flab = parse(text=paste("Total~Fishing~Mortality~(","y^-1",")",sep=""))
#       plot(xrange,yrange,type='n', xlab=Flab, ylab="Total Yield (mt)")
#
#       double.lines(Fmort,predC,bcol="darkgreen",fcol="lightgreen",lwd=5) 
#       points(Fmort,obsC,col= "darkgreen",pch=3,cex=2)
#       points(Fmort,predC,col="darkgreen",pch=16)
#   #   wmaxC = which(obsC==max(obsC,na.rm=TRUE))
#   #   lines(c(0,Fmort[wmaxC]),c(0,obsC[wmaxC]),col="red",lty="longdash")
#       lines(Fyield,yield,col="red",lwd=3,lty="longdash")
#
#       text(x=Fmort[wy5],y=obsC[wy5],labels=floor(dat$t[wy5]),
#             pos=4,offset=0.5,cex=0.8)
#       show.block.number(block,dat$t[1],line=3)
   } #if (plot.prod)

   if (plot.impact)
   {
       d = d + 1
       if (devices[d] > 0)
       {
          s = dev.set(devices[d])
       }
       else
       {
          width = 9.0
          height =9.0
          xpos = xpos + 50
          ypos = ypos + 50
          x11(width=width,height=height,xpos=xpos,ypos=ypos,
               title="Impact")
          devices[d] = dev.cur()
       }
       par(mar=c(4,5,0,4)+0.1)
       ntime = length(dat$t)
       x = dat$t
       Fmort = rowSums(dat[,F.ndx])
       y = matrix(nrow=ntime,ncol=3)
       y[,3] = dat$pop1 + dat$pop2
       y[1,1] = y[1,3]
       y[1,2] = y[1,3]
       for (t in 2:ntime)
       {
          for (p in 1:2)
          {
             FF = Fmort[t]
             if (p == 2)
                FF = 0.0

             y[t,p]  = (K*(r-FF))/((((K*(r-FF))/y[t-1,p])*exp(-(r-FF))) - r*exp(-(r-FF))  + r) # p5
          } 
      }

      legend = c(" F (est)"," F= 0", " B (est)")
      xrange=nice.ts.plot(x,y[,1:2],ylab="Biomass (mt)",legend=legend[1:2])

   #  impact = 100.0*(1.0 - y[,1]/y[,2]) # PNAS
   #  impact = (1.0 - y[,1]/y[,2])
      impact = y[,1]/y[,2]  # WCPFC
   #  print(impact)
      par("new"=TRUE)
      plot(x,impact,type='l',xlim=xrange,ylim=c(0,1),lwd=3,lty="dashed",col="red",
           ann=FALSE,axes=FALSE)
      axis(4,col="red",ylab="p",col.axis="red")
      #mtext("p",side=4,col="red",line=0.1)
      show.block.number(block,dat$t[1],line=3)

   } #if (plot.impact)

   if (plot.Discr)
   {
       d = d + 1
       if (devices[d] > 0)
       {
          s = dev.set(devices[d])
       }
       else
       {
          width = 9.0
          height =11.0
          xpos = xpos + 50
          ypos = ypos + 50
          x11(width=width,height=height,xpos=xpos,ypos=ypos,title="Discriminant")
          devices[d] = dev.cur()
       }

       ntime = length(dat$t)
       gg = c((gear.col+1):(gear.col+ngear))
       print(gg)
       Fmort = rowSums(dat[,gg])
       print(Fmort)
       x = dat$t
       a = -4.0*dat$K #-4.0*r*K
       print("a:")
       print(a)
       b = r - Fmort - T12 - 2.0*q*r*dat$pop1/dat$K  # r - F - T12 - 2.0*q*r*N1/K
       print("b:")
       print(b)
       c = dat$forcing/dat$K #T21/K 
       print("c:")
       print(c)
       Discr = b^2 - 4.0*a*c
       print("Discr:")
       print(Discr)
       plot(x,Discr,type='l')

   } # if (plot.Discr)

   new.devices = devices
   return(new.devices)
}

log.diagnostics=function(file="xssams_program.log",ntime=61,dt=1,ngear=5,
                         plot.Fmort=FALSE,plot.prod=FALSE,plot.impact=FALSE,
                         plot.Discr=FALSE)
{
      
   print(paste("Scanning file",file))
   log = scan(file,what="character")
   res = grep("Residuals:",log)
#   logsdlogF = grep("logsdlogF:",log)   
#   sdlogF = exp(as.numeric(log[logsdlogF+1]))
##  print(sdlogF)
#   logsdlogPop = grep("logsdlogPop:",log)
#   sdlogPop = exp(as.numeric(log[logsdlogPop+1]))
##  print(sdlogPop)
#   logsdlogYield = grep("logsdlogYield:",log)
##  print(length(logsdlogYield))
#   sdlogYield = exp(as.numeric(log[logsdlogYield+1]))
##  print(sdlogYield)
#
#   T12.pos = grep("^T12",log)
##  print(T12.pos)
#   T12 = as.numeric(log[T12.pos+2])
#   wT12 = which(!is.na(T12))
#   T12 = T12[wT12]
##  print(T12)
#
#   K.pos = grep("^K",log)
#   K1 = as.numeric(log[K.pos+2])
#   wK1 = which(!is.na(K1))
#   K = K1[wK1]
##  print(K)
#
#   qProp.pos = grep("^qProp",log)
#   qProp1 = as.numeric(log[qProp.pos+2])
#   wqProp1 = which(!is.na(qProp1))
#   qProp = qProp1[wqProp1]
##  print(qProp)
#
#   r.pos = grep("^r",log)
#   r1 = as.numeric(log[r.pos+2])
#   wr1 = which(!is.na(r1))
#   r = r1[wr1]
##  print(r)
#   sdlogQ.pos = grep("^sdlogQ:",log)
#   sdlogQ = as.numeric(log[sdlogQ.pos+1])
##  print(paste("sdlogQ",sdlogQ))
#  
##  if(1)
##    return(FALSE)
   max.counter = length(res)
   counter = max.counter
#   print(paste(max.counter, "blocks found:"))
##  print(res)
#
#   ncol = (3*ngear+6)
#   diag = matrix(nrow=ntime,ncol=ncol)
#   cnames = vector(length=ncol)

   c = 'n'
   dev.list = vector(mode="numeric",length=6)
   dev.file.names=c("tmp/est_pop","tmp/est_catch","tmp/est_F")
#  while (c != 'q')
   while ( (c != 'q') && (c != 'x') )
   {
      tmp=get.diagnostics(log,ntime=ntime,dt=dt,ngear=ngear,block=counter,mtype="x")
#      fc = res[counter]
#      for (c in 1:ncol)
#      {
#         fc = fc + 1
#      #  print(paste(fc,log[fc])) 
#         cnames[c] = log[fc]
#      }
#      colnames(diag) = cnames
#   
#      for (t in 1:ntime)
#      {
#         for (c in 1:ncol)
#         {
#            fc = fc + 1
#            diag[t,c] = as.numeric(log[fc],ngear)
#         }
#      }
##     print(paste("Block:",counter))
##     print(head(diag))
##     print(tail(diag))
#      print(paste("Displaying block ",counter,sep=""))
#      new.devices = plot.diagnostics(as.data.frame(diag),dt=dt,ngear=ngear,
#                    sdlogPop=sdlogPop[counter], 
#                    sdlogYield=sdlogYield[counter], 
#                    sdlogF=sdlogF[counter],
#                    sdlogQ=sdlogQ[counter],
#                    K=K[counter], r=r[counter],T12=T12[counter],
#                    q=qProp[counter],
#                    plot.Fmort=plot.Fmort,
#                    plot.prod=plot.prod,
#                    plot.impact=plot.impact,
#                    plot.Discr=plot.Discr,
#                    devices=dev.list,block=counter)
      new.devices = plot.diagnostics(tmp$resid,dt=dt,ngear=ngear,
                    sdlogPop=tmp$sdlogPop, 
                    sdlogYield=tmp$sdlogYield, 
                    sdlogF=tmp$sdlogF,
                    sdlogQ=tmp$sdlogQ,
                    K=tmp$K, r=tmp$r,T12=tmp$T12,
                    q=tmp$qProp,
                    plot.Fmort=plot.Fmort,
                    plot.prod=plot.prod,
                    plot.impact=plot.impact,
                    plot.Discr=plot.Discr,
                    devices=dev.list,block=counter)

      dev.list=new.devices
     
      c = readline("next, back, save, quit or exit? (enter n,b,s,q,x):")
      print(paste(c," entered"))
      if (c == 'n')
         counter = counter + 1
      else if (c == 'b')
         counter = counter - 1
      else if (c == 's')
      {
         print(dev.list)
         for (d in 1:length(dev.list))
         {
            if (dev.list[d] > 0)
            {
               tname = paste(dev.file.names[d],"B",counter,sep="")
               print(paste("Saving device",d,"to file",tname))
               dev.set(dev.list[d])
               din = par("din")
               print(din)
               save.png.plot(tname,width=din[1],height=din[2])
            }
            else
               print(paste("Skipping device",d))
          }
      }
      if (counter < 1)
         counter = max.counter
      else if (counter > max.counter)
         counter = 1
   } #while
   graphics.off()
   if (c == 'x')
      q("no")
#  return(counter)
}

#show.block.number=function(block.number,x,line=3)
#{
#   mtext(text=paste("(",block.number,")",sep=""),side=1,line=line,
#          at=c(x,0),cex=0.8)
#}

plot.prod.curve=function(file="xssams.rep")
{
   rep=read.table(file=file,header=TRUE)
   pop.ndx=grep("pop",names(rep)) 
   F.ndx=grep("F",names(rep)) 
   predC.ndx=grep("predC",names(rep)) 
   obsC.ndx=grep("obsC",names(rep)) 
   print(names(rep[,F.ndx]))
   print(names(rep[,obsC.ndx]))

   Fmort = rowSums(exp(rep[,F.ndx]))
   obsC = rowSums(exp(rep[,obsC.ndx]))
   plot(Fmort,obsC)
   lines(c(0,max(Fmort,na.rm=TRUE)),c(0,max(obsC,na.rm=TRUE)))
}

plot.resid.hist=function(dat)
{


}

do.rep.junk=function()
{
   rep=read.table("xssams.rep",header=TRUE)
   head(rep)
   dim(rep)
   rep[,1]
   names(rep)
   names(rep[12:16])
   PredC=exp(rep[12:16])
   rowSums(PredC)->PredCt
   plot(PredCt)
   
   FF=c(7:11)
   head(rep[,FF])
   expF=exp(rep[,FF])
   Ft=rowSums(expF)
   plot(Ft)
   plot(PredCt/Ft)
   plot(Ft,PredCt,type='b')

}

