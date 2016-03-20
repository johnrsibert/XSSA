
# 1 pound = 0.00045359237 metric tons
tplb <- 0.00045359237

yq.time = function(year,quarter)
{
   time = year + (quarter-1)*0.25+0.125
   return(time)
}

make.hdar.dat<-function(file="hdar_calendar.csv",yr1=1952,yr2=2012,ngear=5)
{
   print(paste("Reading ",file))
   dat1 <- read.csv(file,header=TRUE)
#  print(head(dat1))
#  print(tail(dat1))
   print(dim(dat1))

   g1 <- unique(dat1$Method)
   print(g1)
   g2 <- as.vector(c("Tuna HL","Troll","Longline","Bottom/inshore HL","Aku boat","Misc"),mode="character")
   print(g2)
   # tuna handline, troll, longline, bottom fish/inshore handline and aku boat

   qname <- as.vector(c("JAN_MAR","APR_JUN","JUL_SEP","OCT_DEC"))
#  print(qname)

   nrd2 = (yr2-yr2+1)*4
   dat2 <- as.data.frame(matrix(nrow=nrd2,ncol=9))
   colnames(dat2)<-c("year","quarter","time","TunaHL","Troll","Longline","InshoreHL","AkuBoat","Misc")
   print(dim(dat2))

   m <- nrow(dat1)
   for (i in 1:m)
   {
      year <- dat1[i,]$Year
      if ( (year >= yr1) && (year < (yr2+1)) )
      {
      #  print(paste(dat1[i,]$Method,get.index(dat1[i,]$Method,g2,nomatch=6)))
         gi <- get.index(dat1[i,]$Method,g2,nomatch=6)

         qi <- get.index(dat1[i,]$Quarter,qname)

         time = year + (qi-1)*0.25+0.125
         j = (year-yr1)*4+qi
      #  print(paste(year,qi,j,time))

         dat2[j,(gi+3)] <- tplb*dat1[i,]$Lbs
         dat2[j,3] <- time
         dat2[j,2] <- qi
         dat2[j,1] <- year
      }
   }
#  print("dat2:")
#  print(dim(dat2))
#  print(head(dat2))
#  print(tail(dat2))

   if (ngear == 5)
   {
#     makes five gear data
      cc = c(1,3:8)
#  print(cc)
      dat3 <- as.matrix(dat2[,cc])
   }
   else if (ngear == 4)
   {
#     makes four gear data by combining "TunaHL" and "InshoreHL"
      for (j in 4:8)
      {
         w = which(is.na(dat2[,j]))
         if (length(w) > 0)
            dat2[w,j] = 0.0

      }
      dat3 = matrix(nrow=nrow(dat2),ncol=6)
      dat3[,1] = dat2[,1]
      dat3[,2] = dat2[,3]
      dat3[,3] = dat2[,4] + dat2[,7]
   #  dat3[,3] = sum(dat2[,4],dat2[,7],na.rm=TRUE)
      dat3[,4] = dat2[,5]
      dat3[,5] = dat2[,6]
      dat3[,6] = dat2[,8]
      colnames(dat3)=c("year","time","Handline","Troll","Longline","Aku Boat")
   }
   else
   {
      print(paste(ngear,"is an icorrect number of gears; you must specify 5 or 4")) 
      return(NULL)
   }
#  print("dat3:")
#  print(dim(dat3))
#  print(head(dat3))
#  print(tail(dat3))
   
   dat.file = paste("hdar_",yr1,"_",yr2,".dat",sep="")
   print(paste("writing",dat.file))
   write(dat3,file=dat.file,ncolumns=dim(dat3)[1])
   print(paste("finished",dat.file))

   return(as.data.frame(dat3))
}


get.index<-function(pattern,x,nomatch = 0)
{
    y <- grep(pattern,x,fixed=TRUE)
    if (length(y) == 1)
       return(y)
    else
       return(nomatch)
}

plot.hdar<-function(dat)
{
   xrange<-c(1948,2015)
   yrange<-c(0,500) #1300000)

   plot(xrange,yrange,type='n')

   for (j in 4:ncol(dat))
   {
      lines(dat[,3],dat[,j],col=j,lwd=3)

   }
   print(colnames(dat))
   ltext <- c(colnames(dat)[4:9])
   print(ltext)

   legend(xrange[1],yrange[2],ltext,lwd=3,col=c(4:ncol(dat)))

}

###############################
plot.all.ccf<-function(dat)
{
   width <- 6.5
   height <- 9.0
   x11(width=width,height=height)
   old.par <- par(no.readonly = TRUE) 
   par(mar=c(2,2.5,0,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   ng = ncol(dat) - 3
   np = ng*(ng-1)/2
   lm <- layout(matrix(c(1:(np+1)),ncol=2,byrow=TRUE))
   layout.show(lm)
   tmp = dat
   for (j in 4:ncol(dat))
   {
      w = which(is.na(tmp[,j]))
      if (length(w) > 0)
         tmp[w,j] = 0.0
   }


   n = 0
   for (j in 4:(ncol(dat)-1))
   {
      dj = diff(tmp[,j])
      for (k in (j+1):(ncol(dat)))
      {
         n = n+1
         print(paste(j,k,n))
         dk = diff(tmp[,k])
       # ccf(tmp[,j],tmp[,k],main="",lwd=2,ylab="")
         ccf(dj,dk,main="",lwd=2,ylab="")
         abline(v=0.0,col="red",lwd=1,lty="dashed")
         title(paste(colnames(tmp)[j],"vs",colnames(tmp)[k]),line=-1,
               font.main=1)
      }
   }


   save.png.plot("ccf",width=width,height=height)
   par(old.par)
}

plot.all.pacf<-function(dat)
{
   width <- 6.5
   height <- 9.0
   x11(width=width,height=height)
   old.par <- par(no.readonly = TRUE) 
   par(mar=c(2,4,0,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   np = ncol(dat) - 3
   lm <- layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)

   for (j in 4:ncol(dat))
   {
      tmp = dat[,j]
      w = which(is.na(tmp))
      tmp[w] = 0.0
      pacf(tmp,main="",ylab=colnames(dat)[j],lwd=3)
   }
#  plot(c(1,23),c(-0.5,0.5),type='n',ylab="HI LL")
#  mtext("Comming soon: More Longline Data",side=3,line=-3)
   save.png.plot("partial_acf",width=width,height=height)
   par(old.par)

}

plot.all.hdar.gears<-function(dat)
{
   width <- 6.5
   height <- 9.0
   x11(width=width,height=height)
   old.par <- par(no.readonly = TRUE) 
   par(mar=c(2,2.75,0,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   np = ncol(dat) - 3
   lm <- layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)

   for (j in 4:ncol(dat))
   {
      nice.ts.plot(dat[,3],dat[,j],lwd=3,label=colnames(dat)[j],
                   ylab="Catch (mt)")
      wna = which(is.na(dat[,j]))
      zz = vector(length=length(wna))
      points(dat[wna,3],zz,pch='|',col="red")
      dy = diff(dat[,j])
      maxdy=max(abs(dat[,j]))
      yrange=c(-maxdy,maxdy)
      par("new"=TRUE)
      plot(dat[-1,3],dy,type='l',axes=FALSE,ann=FALSE,
           lty="dotted",lwd=2,col="red")
   }

   save.png.plot("hdar_catch_history",width=width,height=height)
   par(old.par)
}

get.SD<-function(Yield)
{
  print("Yield",quote=FALSE)
  print(summary(Yield[,4:8]),quote=FALSE)
  logYield = as.matrix(log(dat[,4:8]))
  print("logYield",quote=FALSE)
  print(summary(logYield),quote=FALSE)
  diff.logYield = diff(logYield)
  print("diff.logYield",quote=FALSE)
  print(summary(diff.logYield),quote=FALSE)
  print("median(diff.logYield)",quote=FALSE)
  print(as.vector(summary(diff.logYield)[3,]),quote=FALSE)
  ngear = ncol(diff.logYield)
# md.diff.logYield = vector(length=ngear)
# for (g in 1:ngear)
#   md.diff.logYield[g] = summary(diff.logYield)[3,g]
# print("log(median(diff.logYield))",quote=FALSE)
# print(md.diff.logYield)
# print(log(md.diff.logYield))
  sd.diff.logYield = vector(length=ngear)
  for (g in 1:ngear)
     sd.diff.logYield[g] = sd(diff.logYield[,g],na.rm=TRUE)
  print("sd.diff.logYield",quote=FALSE)
  print(sd.diff.logYield,quote=FALSE)
  print("log(sd.diff.logYield)",quote=FALSE)
  print(log(sd.diff.logYield),quote=FALSE)
  print("mean(log(sd.diff.logYield))",quote=FALSE)
  print(mean(log(sd.diff.logYield)),quote=FALSE)
}

