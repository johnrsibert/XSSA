source("/home/jsibert/Projects/xssa/scripts/utils.R")

plot.diagnostics=function(dat=NULL,file="diagnostics.dat",ngear,devices)
{
   if (is.null(dat))
     dat = read.table(file=file,header=TRUE)

   if (dat$t[1] == 1)
   {
      dat$t = (1951.875 + dat$t*0.25)
      print(names(dat))
   }
   ncol = ncol(dat)
   dd = c(2:3,6:ncol)
#  print(dd)
   print(names(dat)[dd])
   dat[,dd] = exp(dat[,dd])

#  g2 <- as.vector(c("Tuna HL","Troll","Longline","Bottom/inshore HL","Aku boat","Misc"),mode="character")
   # tuna handline, troll, longline, bottom fish/inshore handline and aku boat


   gear.names = c("TunaHL","Troll","Longline","Bottom/inshore HL","AkuBoat")
   title.line = -1
   lwd = 3

   width = 9.0
   height = 11.0
   d = 1
   if (devices[d] > 0)
   {
      s = dev.set(devices[d])
   }
   else
   {
      x11(width=width,height=height)
      devices[d] = dev.cur()
   }

   old.par = par(no.readonly = TRUE) 
   par(mar=c(3,4.5,0,0)+0.1)
   np = 3
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)

   nice.ts.plot(dat$t,dat$pop1,lwd=lwd)
   title(main="Population 1",line=title.line)

   nice.ts.plot(dat$t,dat$pop2,lwd=lwd)
   title(main="Population 2",line=title.line)

   nice.ts.plot(dat$t,dat$K,bcol="black",fcol="lightgray",lwd=lwd-2)
   double.lines(dat$t,(dat$pop1+dat$pop2),bcol="blue",fcol="lightblue",lwd=lwd)
   par("new"=TRUE)
   plot(dat$t,dat$propL,type='l',ann=FALSE,axes=FALSE,
        col="red4",lwd=2,lty="solid",ylim=c(0,1))
   axis(4,line=0,outer=FALSE,tcl=0.5,labels=FALSE,col="red")
   abline(h=0.9,col="red1",lty="dotdash")
   title(main="Total Population",line=title.line)

   width = 9.0
   height =11.0
   d = 2
   if (devices[d] > 0)
   {
      s = dev.set(devices[d])
   }
   else
   {
     x11(width=width,height=height)
     devices[d] = dev.cur()
   }

   par(mar=c(3,4.5,0,0)+0.1)
   np = ngear
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   for (g in 1:ngear)
   {
      nice.ts.plot(dat$t,dat[,(5+ngear+g)],bcol="darkgreen",fcol="lightgreen",lwd=lwd)
      points(dat$t,dat[,(5+2*ngear+g)],col= "darkgreen",pch=16)
   #  title(main=paste("Catch, gear",g))
      title(main=paste("Catch,",gear.names[g]),line=title.line)
   }

   width = 9.0
   height =11.0
   d = 3
   if (devices[d] > 0)
   {
      s = dev.set(devices[d])
   }
   else
   {
      x11(width=width,height=height)
      devices[d] = dev.cur()
   }
   par(mar=c(3,4.5,0,0)+0.1)
   np = ngear
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   for (g in 1:ngear)
   {
      nice.ts.plot(dat$t,dat[,(g+5)],bcol="orange4",fcol="orange",lwd=lwd)
   #  title(main=paste("F mort, gear",g))
      title(main=paste("F mort,",gear.names[g]),,line=title.line)
   }

   new.devices = devices
   return(new.devices)
}

log.diagnostics=function(file="xssams_program.log",ntime=244,ngear=5)
{
      
   print(paste("Scanning file",file))
   log = scan(file,what="character")
   res = grep("Residuals:",log)
   afters = grep("after",log)

   max.counter = length(res)
   counter = max.counter
   print(paste(max.counter, "blocks found:"))
   print(res)

   ncol = (3*ngear+5)
   diag = matrix(nrow=ntime,ncol=ncol)
   cnames = vector(length=ncol)

   c = 'n'
   dev.list = vector(mode="numeric",length=3)
   dev.file.names=c("est_pop","est_catch","est_F")
   while (c != 'q')
   {
      fc = res[counter]
      for (c in 1:ncol)
      {
         fc = fc + 1
      #  print(paste(fc,log[fc])) 
         cnames[c] = log[fc]
      }
      colnames(diag) = cnames
   
      for (t in 1:ntime)
      {
         for (c in 1:ncol)
         {
            fc = fc + 1
            diag[t,c] = as.numeric(log[fc],ngear)
         }
      }
#     print(paste("Block:",counter))
#     print(head(diag))
#     print(tail(diag))

      print(paste("Displaying block ",counter,sep=""))
      new.devices = plot.diagnostics(as.data.frame(diag),ngear=ngear,devices=dev.list)
      dev.list=new.devices
     

#     title(main=paste(log[afters[counter]+1],"entries"))

      c = readline("next, back, save or quit? (enter n,b,s,or q):")
      print(paste(c," entered"))
      if (c == 'n')
         counter = counter + 1
      else if (c == 'b')
         counter = counter - 1
      else if (c == 's')
      {
         for (d in 1:length(dev.list))
         {
            dev.set(dev.list[d])
            din = par("din")
            save.png.plot(dev.file.names[d],width=din[1],height=din[2])
          }
      }
      if (counter < 1)
         counter = max.counter
      else if (counter > max.counter)
         counter = 1
   } #while
   graphics.off()
#  return(counter)
}



#plot.diagnostics.1=function(dat=NULL,file="diagnostics.dat")
#{
#   if (is.null(dat))
#     dat = read.table(file=file,header=TRUE)
#
#   if (dat$t[1] == 1)
#   {
#      dat$t = (1951.875 + dat$t*0.25)
#   #  print(head(dat))
#   }
##  prop = exp(dat$pop1)/(exp(dat$pop1)+exp(dat$pop2))
##  print(summary(cbind(dat,prop)))
#   dat$pop1 = exp(dat$pop1)
##  w = which(dat$pop1 <= 2)
##  print(w)
##  is.na(dat$pop1[w])
#   dat$pop2 = exp(dat$pop2)
##  print(head(cbind(dat$pop1,dat$pop2))) 
##  dat$predC1=exp(dat$predC1)
##  dat$obsC1=exp(dat$obsC1)
#
#   width = 6.5
#   height = 9.0
#   x11(width=width,height=height)
#   print(paste("Current device",dev.cur()))
#   old.par = par(no.readonly = TRUE) 
#   par(mar=c(3,4.5,0,0)+0.1)
#   np = 3
#   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
#   layout.show(lm)
#   lwd = 3
#
#   nice.ts.plot(dat$t,cbind(dat$pop1,dat$pop2),bcol="blue",fcol="lightblue",lwd=lwd)
#   par("new"=TRUE)
#   plot(dat$t,dat$propL,ylim=c(0,1),type='l',col="red",lwd=lwd,axes=FALSE,ann=FALSE)
#
#   nice.ts.plot(dat$t,exp(dat$F1),bcol="darkgreen",fcol="lightgreen",lwd=lwd)
#
#   nice.ts.plot(dat$t,cbind(dat$predC1,dat$obsC1),bcol="orange4",fcol="orange",lwd=lwd)
#
#   return(dat)
##  return(cbind(dat$pop1,dat$pop2))
#}

