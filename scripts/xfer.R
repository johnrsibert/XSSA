######### mostltly logistic stuff below:
require("mvtnorm")
test.bvnorm=function(sx,sy,rho,n=100)
{
   sigma=matrix(ncol=2,nrow=2)
   sigma[1,1]=sx
   sigma[2,2]=sy
   sigma[1,2]=rho*sx*sy
   sigma[2,1]=sigma[1,2]
   print(sigma)
#  x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma, method="cho
   p = rmvnorm(n=n,mean=c(0,0),sigma=sigma,method="chol")
   print(var(p))
   print(cov(p))
   print(cor(p))
   print(summary(p))
   plot(p)

   return(p)
}

plot.logistic.xfer=function(r=0.5,K=1.0,F1=r/2,T12=0.0,T21=0.0,std=c(0.0,0.0,0.0),q=0.1,dt=0.5,ntime=100)
{

   main = paste("F=",F1,", T12=",T12,", T21=",T21,", S=(",std[1],",",std[2],",",std[3],")",", q=",q,sep="")
   ss = 1.0/dt
   sub = paste("r = ",r,", K= ",K,", dt = ",dt,", (ss = ",ss,")",sep="")
   print(main)
   print(sub)
   Fmsy = r/2
   Nmsy = K*(1-Fmsy/r)
   MSY = Fmsy*Nmsy

   res = logistic.xfer(r,K,F1,T12,T21,std,dt,q,ntime)

   xrange = range(res$T)
   yrange = c(0.0,max(res$NN,K))
   lwd  = 5
   lwdab  = lwd-2 
   plot(res$T,res$NN,xlim=xrange,ylim=yrange,type='b',las=1,
         main=main,sub=sub,xlab="Time",ylab="N or C")
   lines(res$T,res$CNN,col="red",lwd=lwd)
   lines(res$T,res$N11,col="blue",lty="dashed",lwd=lwd)
   lines(res$T,res$N21,col="blue",lty="dotted",lwd=lwd)
   lines(res$T,(res$N11+res$N21),col="blue",lwd=lwd)
   lines(res$T,res$N11/(res$N11+res$N21),col="green",lwd=lwd)
   abline(h=K,lty="dotdash",lwd=lwdab)
   abline(h=Nmsy,lty="dotdash",lwd=lwdab)
   abline(h=MSY,lty="dotdash",col="red",lwd=lwdab)
   abline(h=0.9,lty="dotdash",col="green",lwd=lwdab)

   file = paste("r",r,"F",F1,"T12",T12,"T21",T21,"S",std[1],std[2],std[3],"q",q,sep="")
   print(file)
   file = gsub(".","",file,fixed=TRUE)
   print(file)
   save.png.plot(file)

   return(res)
}

draw.examples=function()
{
   plot.logistic.xfer()->tmp
   x11()
   plot.logistic.xfer(T21=0.1)->tmp
   x11()
   plot.logistic.xfer(T21=0.01)->tmp
   x11()
   plot.logistic.xfer(T21=0.001)->tmp
   x11()
   plot.logistic.xfer(T12=0.1)->tmp
   x11()
   plot.logistic.xfer(T12=0.01)->tmp
   x11()
   plot.logistic.xfer(T12=0.001)->tmp
   x11()
   plot.logistic.xfer(T21=0.001,T12=0.01)->tmp
   x11()
   plot.logistic.xfer(T21=0.001,T12=0.01,dt=0.1,std=c(0.025,0.025,-0.5))->tmp
}

logistic.xfer=function(r,K,F1,T12,T21,cov,dt,q,ntime)
{
   N11=vector(length=ntime)
   N11[1] = K
   N21=vector(length=ntime)
   N21[1] = 0.0
   NN=vector(length=ntime)
   NN[1] = N11[1]+N21[1]

   ss = 1.0/dt
   sigma=matrix(ncol=2,nrow=2)
   sigma[1,1]=cov[1]
   sigma[2,2]=cov[2]
   sigma[1,2]=cov[3]*cov[1]*cov[2]
   sigma[2,1]=sigma[1,2]
   if (sum(abs(cov)) > 0.0)
      std = rmvnorm(n=ntime,mean=c(0,0),sigma=sigma,method="chol")
   else
      std = matrix(0.0,ncol=2,nrow=ntime)
 
#  print(paste(p[1],T12,T21))
   for (i in 2:ntime)
   {
      prevN11 = N11[i-1]
      prevN21 = N21[i-1]
      prevNN  = NN[i-1]
   #  print(paste(p[i],T12,T21))
      for (k in 1:ss)
      {
         nextNN = prevNN + dt*r*prevNN*(1.0 - prevNN/K) - dt*prevNN*(F1+T12) + dt*T21
         prevNN = nextNN

         nextN11 = prevN11 + dt*r*prevN11*(1.0 - prevN11/K) - dt*prevN11*(F1+T12) - 2.0*(1.0-q)*r*dt*prevN11*prevN21/K
         prevN11 = nextN11
         nextN21 = prevN21 + dt*r*prevN21*(1.0 - prevN21/K) - dt*prevN21*(F1+T12) - 2.0*q*r*dt*prevN11*prevN21/K + dt*T21
         prevN21 = nextN21

      }
      NN[i] = nextNN
      N11[i] = nextN11*exp(std[i,1])
      N21[i] = nextN21*exp(std[i,2])
   }

#  CNN<-F1*NN
   CNN<-F1*(N11+N21)

   ret = list(T=seq(1:ntime),NN=NN, N11=N11, N21=N21, CNN=CNN)
   return(ret)
}

plot.logistic=function(r=0.5,K=1.0,dt=1)
{
   ntime = 100
   N=matrix(nrow=ntime+1,ncol=2)
   tt=seq(0:ntime)
   N0=vector(length=2)
   N0[1]= 0.75 * K #1e-8
   N0[2] = 1.25*K #2.0*K
   N[1,1] = N0[1]
   N[1,2] = N0[2]
   ss = 1.0/dt
   dtr = r*dt
   title = paste("r = ",r,", K= ",K,", dt = ",dt,", (ss = ",ss,")",sep="")

   #exact
   for (i in 2:ntime)
   {
      exprt = exp(r*tt[i])
      for(n in 1:2)
      {
        N[i,n] = N0[n]*exprt/(1.0 - (N0[n]/K)*(1-exprt))
      }
   }
   plot(tt,N[,1],type='b',ylim=c(N0[1],N0[2]),main=title,ylab="N",xlab="t",bty="l")
   lines(tt,N[,2],type='b')

   # explicit (1)
   for (i in 2:ntime)
   {
      for(n in 1:2)
      {
         prevN = N[i-1,n]
         for (k in 1:ss)
         {
            nextN = prevN + dtr*prevN*(1.0 - prevN/K)
            prevN = nextN
         }
         N[i,n] = nextN
      }
   }  
   lines(tt,N[,1],type='l',col="red",lwd=2)
   lines(tt,N[,2],type='l',col="red",lwd=2)


   # partially implicit (2a)
   for (i in 2:ntime)
   {
      for(n in 1:2)
      {
         prevN = N[i-1,n]
         for (k in 1:ss)
         {
            nextN = prevN/(1.0 - dtr*(1.0 - prevN/K))
            prevN = nextN
         }
         N[i,n] = nextN
      }
   }  
   lines(tt,N[,1],type='l',col="blue",lwd=2)
   lines(tt,N[,2],type='l',col="blue",lwd=2)

   # partially implicit (2b)
   for (i in 2:ntime)
   {
      for(n in 1:2)
      {
         prevN = N[i-1,n]
         for (k in 1:ss)
         {
            nextN = prevN*(1.0+dtr)/(1.0 + dtr*prevN/K)
            prevN = nextN
         }
         N[i,n] = nextN
      }
   }  
   lines(tt,N[,1],type='l',col="green",lwd=2)
   lines(tt,N[,2],type='l',col="green",lwd=2)

   # fully implicit (3)
   a = r/K
   b = 1-r
   for (i in 2:ntime)
   {
      for(n in 1:2)
      {
         prevN = N[i-1,n]
         for (k in 1:ss)
         {
            c = -prevN
            rad = sqrt(b*b - 4.0*a*c)
            root1 = (-b + rad)/(2.0*a)
        #   root2 = (-b - rad)/(2.0*a)
            prevN = root1
         }
      #  print(paste(root1,root2,sep=", "))
         N[i,n] = root1
      }
   }  
   lines(tt,N[,1],type='l',col="orange",lwd=2)#,lty="dashed")
   lines(tt,N[,2],type='l',col="orange",lwd=2)#,lty="dashed")

   file = paste("r",r,"K",K,"dt",dt,"ss",ss,sep="")
   print(file)
   file = sub(".","",file,fixed=TRUE)
   print(file)
   save.png.plot(file)
}

save.png.plot<-function(root,width=6.5,height=4.5)
{
  graphics.root <-paste("./Reports/graphics/",root,sep="")
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

