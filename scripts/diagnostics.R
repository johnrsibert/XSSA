source("/home/jsibert/Projects/xssa/scripts/utils.R")

plot.diagnostics=function(dat=NULL,file="diagnostics.dat",dt,ngear,
                 sdlogPop, sdlogYield, sdlogF,plot.Fmort,
                 devices,block)
{
   print(paste("Block: ",block))
   if (is.null(dat))
     dat = read.table(file=file,header=TRUE)

   start.year = 1952
   if (dat$t[1] == 1)
   {
      dat$t = (start.year-0.4*dt +  dat$t*dt)
   #  print(names(dat))
   }
   ncol = ncol(dat)
#  print("Names of all variables:")
#  print(names(dat))
   gear.col = 6
   dd = c(2:3,(gear.col+1):ncol)
#  print("Names of log transformed variables:")
#  print(names(dat)[dd])
   dat[,dd] = exp(dat[,dd])

   gear.names = c("TunaHL","Troll","Longline","Bottom/inshore HL","AkuBoat")
   title.line = -1
   lwd = 3

   old.NN.plot = FALSE
   if(old.NN.plot)
   {
      width = 9.0
      height =11.0
      x11(width=width,height=height)
     
      old.par = par(no.readonly = TRUE) 
      par(mar=c(3,4.5,0,0)+0.1)
      np = 3
      lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
      layout.show(lm)
    
      nice.ts.plot(dat$t,dat$pop1,lwd=lwd)
      title(main="Population 1",line=title.line)
    
      nice.ts.plot(dat$t,dat$pop2,lwd=lwd)
      title(main="Population 2",line=title.line)
   
      nice.ts.plot(dat$t,dat$K,bcol="black",fcol="lightgray",lwd=lwd)
      double.lines(dat$t,(dat$pop1+dat$pop2),bcol="blue",fcol="lightblue",lwd=lwd)
      title(main="Total Population",line=title.line)
    
      par("new"=TRUE)
      plot(dat$t,dat$forcing,type='n',ann=FALSE,axes=FALSE,
           ylim=c(0.0,max(dat$forcing)))
      double.lines(dat$t,dat$forcing,bcol="purple4",fcol="purple1",lwd=lwd)
      axis(4,line=-1,outer=FALSE,labels=FALSE,tcl=0.5,col="purple1")
       
      par("new"=TRUE)
      plot(dat$t,dat$propL,type='l',ann=FALSE,axes=FALSE,
           col="red4",lwd=2,lty="solid",ylim=c(0,1))
      axis(4,line=0,outer=FALSE,tcl=0.5,labels=FALSE,col="red")
      abline(h=0.9,col="red4",lty="dotdash",lwd=2)
   } #if(old.NN.plot)
   else
   {
      d = 1
      if (devices[d] > 0)
      {
         s = dev.set(devices[d])
      }
      else
      {
         width = 9.0
         height = 9.0
         x11(width=width,height=height)
         devices[d] = dev.cur()
      }

      x = dat$t
      ntime=length(x)
      y = matrix(nrow=length(x),ncol=3)
      y[,1] = dat$pop1
      y[,2] = dat$pop2
      y[,3] = dat$pop1 + dat$pop2
      legend = c(" N1"," N2"," N1+N2")
      xrange=nice.ts.plot(x,y,legend=legend,lwd=5,ylab="Biomass (mt)")
      lines(x,dat$K,lwd=2,lty="dotdash",col="blue",xlim=xrange)
      sdy = exp(log(y[,3])+2.0*sdlogPop)
      lines(dat$t,sdy,col="blue",lty="dotted")
      sdy = exp(log(y[,3])-2.0*sdlogPop)
      lines(dat$t,sdy,col="blue",lty="dotted")

      par("new"=TRUE)
      plot(x,dat$propL,lwd=3,type='l',col="red",ylim=c(0,1),
           ann=FALSE,axes=FALSE,xlim=xrange)
      text(x[ntime],dat$propL[ntime]," p",adj=c(0,0),col="red")
      abline(h=0.9,lwd=2,lty="dotdash",col="red")
      axis(4,col="red",ylab="p",col.axis="red")
      mtext("p",side=4,col="red",line=0.1)

      par("new"=TRUE)
      plot(x,dat$forcing,lwd=3,type='l',col="purple", 
           ylim=c(0.0,max(dat$forcing)), ann=FALSE,axes=FALSE,xlim=xrange)
      text(x[ntime],dat$forcing[ntime]," T21",adj=c(0,0),col="purple")
      axis(4,line=-2,col="purple",col.axis="purple")
   #  axis(2,col="purple",col.axis="purple",pos=xrange[2])
      mtext("T21",side=4,col="purple",line=-2.0)

      title(main=paste("Block", block),line=title.line)

   }

   d = 2
   if (devices[d] > 0)
   {
      s = dev.set(devices[d])
   }
   else
   {
     width = 9.0
     height =11.0
     x11(width=width,height=height)
     devices[d] = dev.cur()
   }

   par(mar=c(3,4.5,0,0)+0.1)
   np = ngear
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   for (g in 1:ngear)
   {
      nice.ts.plot(dat$t,dat[,(gear.col+ngear+g)],bcol="darkgreen",fcol="lightgreen",lwd=lwd,ylab="Catch (mt)")
      points(dat$t,dat[,(gear.col+2*ngear+g)],col= "darkgreen",pch=16)
      if (g == 1)
         title(main=paste(gear.names[g]," (",block,")",sep=""),line=title.line)
      else
         title(main=gear.names[g],line=title.line)
      sdy = exp(log(dat[,(gear.col+ngear+g)])+2.0*sdlogYield)
      lines(dat$t,sdy,col="darkgreen",lty="dotted")
      sdy = exp(log(dat[,(gear.col+ngear+g)])-2.0*sdlogYield)
      lines(dat$t,sdy,col="darkgreen",lty="dotted")
   }

#  plot.Fmort = TRUE
   if (plot.Fmort)
   {
       d = 3
       if (devices[d] > 0)
       {
          s = dev.set(devices[d])
       }
       else
       {
          width = 9.0
          height =11.0
          x11(width=width,height=height)
          devices[d] = dev.cur()
       }
       par(mar=c(3,4.5,0,0)+0.1)
       np = ngear
       lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
       layout.show(lm)
       for (g in 1:ngear)
       {
          nice.ts.plot(dat$t,dat[,(gear.col+g)],bcol="orange4",fcol="orange",lwd=lwd)
          if (g == 1)
             title(main=paste(gear.names[g]," (",block,")",sep=""),line=title.line)
          else
             title(main=gear.names[g],line=title.line)
          sdy = exp(log(dat[,(gear.col+g)])+2.0*sdlogF)
          lines(dat$t,sdy,col="orange4",lty="dotted")
          sdy = exp(log(dat[,(gear.col+g)])-2.0*sdlogF)
          lines(dat$t,sdy,col="orange4",lty="dotted")
       }
   } #if (plot.Fmort)

   new.devices = devices
   return(new.devices)
}

log.diagnostics=function(file="xssams_program.log",ntime=244,dt=0.25,ngear=5,plot.Fmort=FALSE)
{
      
   print(paste("Scanning file",file))
   log = scan(file,what="character")
   res = grep("Residuals:",log)
   logsdlogF = grep("logsdlogF:",log)   
   sdlogF = exp(as.numeric(log[logsdlogF+1]))
#  print(sdlogF)
   logsdlogPop = grep("logsdlogPop:",log)
   sdlogPop = exp(as.numeric(log[logsdlogPop+1]))
#  print(sdlogPop)
   logsdlogYield = grep("logsdlogYield:",log)
   sdlogYield = exp(as.numeric(log[logsdlogYield+1]))
#  print(sdlogYield)
   max.counter = length(res)
   counter = max.counter
   print(paste(max.counter, "blocks found:"))
#  print(res)

   ncol = (3*ngear+6)
   diag = matrix(nrow=ntime,ncol=ncol)
   cnames = vector(length=ncol)

   c = 'n'
   dev.list = vector(mode="numeric",length=3)
   dev.file.names=c("est_pop","est_catch","est_F")
#  while (c != 'q')
   while ( (c != 'q') && (c != 'x') )
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
      new.devices = plot.diagnostics(as.data.frame(diag),dt=dt,ngear=ngear,
                    sdlogPop=sdlogPop[counter], 
                    sdlogYield=sdlogYield[counter], 
                    sdlogF=sdlogF[counter],
                    plot.Fmort=plot.Fmort,
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

plot.resid.hist=function(dat)
{


}

