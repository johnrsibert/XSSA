gear.names = c("TunaHL","Troll","Longline","Bottom/inshore HL","AkuBoat")
sgn = c("THL","Troll","LL","BHL","Aku")
 
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

   nN = vector(length=2)
   nN[1] = exp(log(pN1) + dLN1 * dt)
   nN[2] = exp(log(pN2) + dLN2 * dt)
#  print(paste("n:",nN[1],nN[2]))
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
   image(iN1,iN2,prop,zlim=c(0,1),col=heat.colors(5),add=TRUE)

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

qcomp=function()
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

plot.NN.ts=function(pop,K)
{
   ntime = nrow(pop)

   lprop = pop[,1]/pop[,3]

   width = 9.0
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,4,0,4)+0.1,las=1)

   x = c(1:ntime)
   nice.ts.plot(x,pop,legend=colnames(pop),xlab="t",ylab="N")
   abline(h=K,lwd=2,lty="dotdash",col="blue")
   text(x[ntime],K,"K",adj=c(0,0),col="blue")

   par("new"=TRUE)
   plot(x,lprop,lwd=3,type='l',col="red",ylim=c(0,1),
        ann=FALSE,axes=FALSE)
   abline(h=0.9,lwd=2,lty="dotdash",col="red")
   axis(4,col="red",ylab="p",col.axis="red")
   abline(v=par("usr")[2],lwd=2,col="red")
   mtext("p",side=4,col="red",line=3)
#  title("Log")
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
  
   plot.NN.ts(pop,K)

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
   plot.NN.ts(pop,K)
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


plot.propL.prior=function(propL=0.9)
{
   p = seq(0.01,0.99,0.01)
   p = c(0.0001,p,0.9999)
   wp = which(p == propL)
   
   plot(c(0.0,1.0),c(0.0,1.0),type='n',xlab="x",ylab="p(L(x))")
   abline(v=propL,col="red",lty="dotdash")
   for (sd in c(0.6,0.7,0.8,0.9,0.99,0.999,0.9999))
#  sd = 0.9
   {
      pp = dnorm(logit(p),mean=logit(propL),sd=logit(sd))
      double.lines(p,pp,lwd=5,fcol="lightblue",bcol="blue")
      text(propL,pp[wp],sd,col="blue")
   }

   save.png.plot("propL_prior",width=6.5,height=6.5)
}

# from MHI_sim.R
compute.F<-function(yr1=1952,yr2=2012,cfile="../HDAR/hdar_1952_2012.dat",plot=TRUE)
{
   eps.na = 1e-8
   obs.catch = read.table(file=cfile)
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

xssams.sim=function(r=0.3, K=200000, q=0.54, T12=0.01, T21=0.002,
                  dt = 1.0, s=c(0.0,0.0,0.0), fr=2, p=0.9)
{
   F.matrix = compute.F(plot=FALSE) 
   sum.F = colSums(F.matrix)
   ntime = ncol(F.matrix)

   region.biomass = as.matrix(read.table("../run/total_biomass.dat"))
   immigrant.biomass = T21*region.biomass[fr,]
   pop = matrix(ncol=3,nrow=ntime)
   colnames(pop)=c("N1","N2","N1+N2")
   pop[1,1] = p*K
   pop[1,2] = (1-p)*K
   pop[1,3] = pop[1,1]+pop[1,2]
   for (t in 2:ntime)
   {
      N1 = pop[t-1,1]
      N2 = pop[t-1,2]
      
      tN = step(N1,N2,r,K,sum.F[t],q,T12,immigrant.biomass[t],dt,s)
      pop[t,1:2] = tN
      pop[t,3] = pop[t,1]+pop[t,2]
   }
  
   plot.NN.ts(pop,K)

   plot.sim.catch.ts(pop[,3],F.matrix)
}

plot.sim.catch.ts=function(t.pop,F)
{
   ntime = ncol(F)
   ngear = nrow(F)

   C.by.gear=matrix(nrow=ngear,ncol=ntime)

   for (g in 1:ngear)
   {
      for (t in 1:ntime)
         C.by.gear[g,t] = F[g,t]*t.pop[t]
   }

   total.C = colSums(C.by.gear)
   x = c(1:ntime)
   width = 9.0
   height = 4.5
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(4,4,0,4)+0.1,las=1)
   nice.ts.plot(x,total.C,xlab="t",ylab="C",
                bcol="darkgreen",fcol="lightgreen",lwd=3)
#  points(dat$t,dat[,(5+2*ngear+g)],col= "darkgreen",pch=16)
   par(old.par)
}
