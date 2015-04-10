
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
   print(paste("n:",nN[1],nN[2]))
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
   plot(c(0.0,1.1),c(0.0,1.15),type='n',xlab="N1",ylab="N2")
   image(iN1,iN2,prop,zlim=c(0,1),col=heat.colors(20),add=TRUE)
   pl = 0.9
   p = pN1*(1.0/pl - 1.0)
   double.lines(pN1,p,lwd=15,bcol="darkgreen",fcol="lightgreen")
   lines(c(0.0,K),c(K,0.0), col="red",lwd=3)
   title(main=paste("r=",r,", K=", K, ", F=", FF, ", q=", q,
                    ", T12=", T12, ",T21=", T21,sep=""),
          line=-1,font.main=1,cex=0.8)
   
#  mat = cbind(u,v,df)
#  print(tail(mat))
   vectorField(u,v,xpos=df[,1],ypos=df[,2],vecspec="lonlat",
               headspan=0.05,scale=5,col="blue")
   }
}

qcomp=function(qq=c(0.25,0.5,0.75))
{
   for(k in 1:length(qq))
   {
      x11(width=9,height=9)
      NNphase(q=qq[k])

   }

}
plotN1N2=function(ntime=100,r=0.3, K=1.0, F = 0.007, q=0.54, 
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
      tN = step(N1,N2,r,K,F,q,T12,T21,dt,s)
      pop[t,1:2] = tN
      pop[t,3] = pop[t,1]+pop[t,2]
   #  print(pop[t,])
   }
   lprop = pop[,1]/pop[,3]

   x = c(1:ntime)
   nice.ts.plot(x,pop,legend=colnames(pop))
   par("new"=TRUE)
   plot(x,lprop,lwd=3,type='l',col="red",ylim=c(0,1),
        ann=FALSE,axes=FALSE)
   abline(h=0.9,lty="dotdash",col="red")
   title("Log")
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
   lprop = pop[,1]/pop[,3]

   x = c(1:ntime)
   nice.ts.plot(x,pop,legend=colnames(pop))
   par("new"=TRUE)
   plot(x,lprop,lwd=3,type='l',col="red",ylim=c(0,1),
        ann=FALSE,axes=FALSE)
   abline(h=0.9,lty="dotdash",col="red")
   title("Arithmetic")
}

a.log.disp=function(ntime=100,r=0.3, K=1.0, F = 0.007, q=0.54, 
                  T12=0.01, T21=0.002, p=0.9, 
                  dt = 1.0, s=c(0.0,0.0,0.0))
{

   plotN1N2(ntime, r, K, F, q, T12, T21, p, dt, s)
   x11()
   plotN1N2a(ntime, r, K, F, q, T12, T21, p, dt, s)

}

