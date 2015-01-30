
# 1 pound = 0.00045359237 metric tons
tplb = 0.00045359237

get.index=function(pattern,x,nomatch = 0)
{
    y = grep(pattern,x,fixed=TRUE)
    if (length(y) == 1)
       return(as.numeric(y))
    else
       return(nomatch)
}

make.noaa.dat=function(file="noaa_quarterly.csv",yr1=1995,yr2=2013)
{
   print(paste("Reading ",file))
   dat1 = read.csv(file,header=TRUE)
#  print(head(dat1))
#  print(tail(dat1))
   print(dim(dat1))

   g2 = as.vector(c("Deep Sets","Shallow Sets"))

   nrd2 = (yr2-yr2+1)*4
   dat2 = as.data.frame(matrix(nrow=nrd2,ncol=5))
   colnames(dat2) = c("year","quarter","time","DSLL","SSLL")

   m = nrow(dat1)
   for (i in 1:m)
   {
      year = dat1[i,]$HAULYR
   #  print(paste(i,year))
      if ( (year >= yr1) && (year < (yr2+1)) )
      {
         gi = get.index(dat1[i,]$SET_CATEGORY,g2,nomatch=3)
         qi = dat1[i,]$HAULQTR
#        print(qi)
         time = year + (qi-1)*0.25+0.125
#        print(tt)
         j = (year-yr1)*4+qi
#        print(paste(i,year,qi,j,time))

         dat2[j,(gi+3)] = tplb*dat1[i,]$YFT_LBS
    #    dat2[j,(gi+3)] = dat1[i,]$YFT_LBS
         dat2[j,3] = time
         dat2[j,2] = qi
         dat2[j,1] = year
      }

   }


   return(dat2)
}

LL.join=function(hdar=NULL,noaa=NULL,yr1=1952,yr2=2012)
{
   if(is.null(hdar))
      hdar = make.hdar.dat(file="../HDAR/hdar_calendar.csv",yr1=1949,yr2=2014)
   print(head(hdar))

   if(is.null(noaa))
      noaa = make.noaa.dat(file="../NOAA/noaa_quarterly.csv",yr1=1995,yr2=2013)
   print(head(noaa))

   mm = match(hdar$time,noaa$time)
#  print(paste(nrow(hdar),nrow(noaa),length(mm)))

   ntime = nrow(hdar)
   LL = matrix(nrow=ntime,ncol=6)
   colnames(LL)=c("time","F","D","S","DS","Mean")

   for (n in 1:ntime)
   {
      LL[n,1] = hdar[n,3]
      LL[n,2] = hdar[n,6]
      if (!is.na(mm[n]))
      {
         ni = mm[n]
         LL[n,3] = noaa[ni,4]
         LL[n,4] = noaa[ni,5]
         if (is.na(LL[n,4]))
            LL[n,5] = LL[n,3]
         else if (is.na(LL[n,3]))
            LL[n,5] = LL[n,4]
         else
            LL[n,5] = LL[n,3] + LL[n,4]
      }
   } 
   LL = as.data.frame(LL)
   #rowMeans(cbind(tmp$F,tmp$DS),na.rm=TRUE)
   LL$Mean = rowMeans(cbind(LL$F,LL$DS),na.rm=TRUE)

   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   old.par = par(no.readonly = TRUE) 
   par(mar=c(2,4,0,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   lm <- layout(matrix(c(1:2),nrow=2,byrow=TRUE))
   layout.show(lm)


   plot(LL$time,LL$F,las=1,type='l',lwd=5,col="blue",
              ylab="Yellofin Catch (mt)",xlab="")
   lines(LL$time,LL$F,lwd=3,col="lightblue")
   lines(LL$time,LL$S,col="red",lwd=1,lty="dotted")
   lines(LL$time,LL$D,col="orange",lwd=1)
   lines(LL$time,LL$DS,col="blue",lwd=3)
#  lines(LL$time,LL$Mean,col="green",lwd=1)
   rug(LL$time)

   legend("topleft",col=c("lightblue","blue","red","orange"),
                    lwd=c(5,3,1,1),bty='n',
        legend=c("HDAR Longline","NOAA Total Longline",
                  "NOAA Shallow Set Longline","NOAA Deep Set Longline"))

   dy = diff(LL$Mean)
   plot(LL$time,LL$Mean,las=1,type='l',lwd=5,col="blue",
              ylab="Yellofin Catch (mt)",xlab="")
   lines(LL$time,LL$Mean,lwd=3,col="lightblue")
   par("new"=TRUE)
   plot(LL$time[-1],dy,type='l',lwd=1,lty="dotted",col="black",axes=FALSE,ylab="")
   legend("topleft",col=c("lightblue","black"),
                    lwd=c(5,1),bty='n',lty=c("solid","dotted"),
        legend=c("Average of HDAR and NOAA Total",
                  "First differences"))

   save.png.plot("hdar_noaa_LL_ts",width=width,height=height)

   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   old.par <- par(no.readonly = TRUE) 
   par(mar=c(2,4,0,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   np = 5
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)

   for (j in 2:ncol(LL))
   {
      tmp = LL[,j]
      w = which(is.na(tmp))
      tmp[w] = 0.0
      pacf(tmp,main="",ylab=colnames(LL)[j],lwd=3)
   }
   save.png.plot("LL_partial_acf",width=width,height=height)

   nqd = (yr2-yr1+1)*4
   dat = as.data.frame(matrix(nrow=nqd,ncol=5))
   yy = vector(length=nqd)
   colnames(dat) = c("OffshoreHL","Troll","Longline","InshoreHL","AkuBoat")
   for (n in 1:ntime)
   {
      year = hdar[n,]$year
      if ( (year >= yr1) && (year < (yr2+1)) )
      {
         qi = hdar[n,]$quarter
         j = (year-yr1)*4+qi
         dat[j,1] = hdar[n,4]#$OffshoreHL
         dat[j,2] = hdar[n,5] #$Troll
         dat[j,3] = LL[n,]$Mean # Longline
         dat[j,4] = hdar[n,7]#$InshoreHL
         dat[j,5] = hdar[n,8]#$AkuBoat
         yy[j] = LL$time[n]
      }
   }
   colnames(dat) = c("Tuna HL","Troll","Longline","Bottom/inshore HL","Aku boat")

   dat.file = paste("five_gears_",yr1,"_",yr2,".dat",sep="")
   print(paste("writing",dat.file))
   write(as.matrix(dat),file=dat.file,ncolumns=dim(dat)[1])
   print(paste("finished",dat.file))

   width = 6.5
   height <- 9.0
   x11(width=width,height=height)
   par(mar=c(2.5,4,0,0)+0.1)
#  par(mar=c(5  ,4,4,2)+0.1)
   np = ncol(dat)
   lm <- layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)

   for (j in 1:ncol(dat))
   {
      nice.ts.plot(yy,dat[,j],lwd=3,label=colnames(dat)[j],
                   ylab="Catch (mt)")
   }

   save.png.plot(paste(ncol(dat),"_gear_catch_history",sep=""),width=width,height=height)




   par(old.par)
   return(cbind(yy,dat))
}

make.catch.diffs=function(file="five_gears.dat")
{
   dat = t(read.table(file))
   colnames(dat) = c("OffshoreHL","Troll","Longline","InshoreHL","AkuBoat")
   print(head(dat))

   ZeroCatch = 1.0 #1.0e-8
   logZeroCatch = log(ZeroCatch)
   ngear = ncol(dat)

   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   old.par <- par(no.readonly = TRUE) 
   par(mar=c(2,4,1,0)+0.1)
#  par(mar=c(5,4,4,2)+0.1)
   np = ngear
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   breaks = seq(-4.25,4.25,0.5)

   for (g in 1:ngear)
   {
      ts = dat[,g]
   #  w = which (ts <= 0.0)
   #  ts[w] = ZeroCatch
      ts = log(ts+ZeroCatch)
      dts = diff(ts)
      sd = sd(dts)
      mean = mean(dts)
      print(paste(g,mean,sd))
      hist(dts,breaks=breaks,main=colnames(dat)[g],freq=FALSE,las=1)
      x = seq(-4*sd,4*sd,0.1)
      lines(x,dnorm(x,mean=mean,sd=sd),col="blue")
      lines(x,dtdist(x,mean=mean,sd=1),col="red")
   #  lines(x,dcauchy(x,location=mean,scale=1),col="green" ,lty="dashed")
   }
   save.png.plot("first_difference_histograms",width=width,height=height)
}
#dnorm(x, mean = 0, sd = 1, log = FALSE)
dtdist=function(x, mean = 0, sd = 1)
{
  xx = (x-mean)/sd
# print(paste(xx,xx*xx,(1.0+xx*xx)))
  dt = 1.0/(pi*(1.0+xx*xx))
  return(dt)
}



