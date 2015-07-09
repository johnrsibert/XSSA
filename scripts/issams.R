#require("mvtnorm")
gear.names = c("TunaHL","Troll","Longline","Bottom/inshore HL","AkuBoat")
sgn = c("THL","Troll","LL","BHL","Aku")
have.xssams.R = FALSE
have.xssams.R = TRUE

plot.error=function(x,y,sd,bcol,fcol,mult=2)
{
   sdyu = exp(log(y)+mult*sd)
   sdyl = exp(log(y)-mult*sd)
   polygon(c(x,rev(x)),c(sdyl,rev(sdyu)),
              border=bcol,lty="dashed",lwd=1,col=fcol)

}


plot.diagnostics=function(dat=NULL,file="diagnostics.dat",dt,ngear,
                 sdlogPop, sdlogYield, sdlogF, sdlogQ, K, r,
                 plot.Fmort,plot.prod,
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

   wy5 = which((floor(dat$t%%5))==0)
   ncol = ncol(dat)
#  print("Names of all variables:")
#  print(names(dat))
#  print(head(dat))
#  print(tail(dat))
   gear.col = 4
   dd = c(2,(gear.col+1):ncol)
#  print("Names of log transformed variables:")
#  print(names(dat)[dd])
   dat[,dd] = exp(dat[,dd])
#  print(head(dat))
#  print(tail(dat))

   gear.names = c("TunaHL","Troll","Longline","Bottom/inshore HL","AkuBoat")
   title.line = -1
   lwd = 3
   #sd.lwd = 3
   #sd.lty = "dotted"
   #old.par = par(no.readonly = TRUE) 

   F.ndx=grep("F",names(dat)) 
   predC.ndx=grep("predC",names(dat)) 
   obsC.ndx=grep("obsC",names(dat)) 

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

   par(mar=c(4,5,0,4)+0.1)
   x = dat$t
   ntime = length(dat$t)

   options(scipen=6)
   xrange=nice.ts.plot(dat$t,dat$pop,legend="N",lwd=5,ylab="Biomass (mt)")
#  plot.error(dat$t,dat$pop,sdlogPop,bcol="blue",fcol="lightblue")
   lines(dat$t,dat$pop,col="blue",lwd=5)

   lines(x,dat$K,lwd=2,lty="dotdash",col="blue",xlim=xrange)
   text(x[ntime],dat$K[ntime]," K",adj=c(0,0.5),col="blue")

   par("new"=TRUE)
   nice.ts.plot(dat$t,dat$forcing,legend="I",lwd=5,ylab="")
   plot.error(x,dat$forcing,sdlogQ,bcol="purple4",fcol="purple1")
   lines(x,dat$forcing,col="purple4",lwd=5) 
   text(x[ntime],dat$forcing[ntime],adj=c(0,0),col="purple4")
   axis(4,col="purple4",col.axis="purple4")
   mtext("Index", side=4,col="purple4",line=0.1)
   show.block.number(block,dat$t[1])

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

   par(mar=c(3,4.5,0,0)+0.1)
   np = ngear
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   for (g in 1:ngear)
   {
      nice.ts.plot(dat$t,dat[,(gear.col+ngear+g)],bcol="darkgreen",fcol="lightgreen",lwd=lwd,ylab="Catch (mt)")
      plot.error(dat$t,dat[,(gear.col+ngear+g)],sdlogYield,
                 bcol="darkgreen",fcol="lightgreen")
      lines(dat$t,dat[,(gear.col+ngear+g)],col="darkgreen",lwd=lwd+2)
      points(dat$t,dat[,(gear.col+2*ngear+g)],col= "darkgreen",pch=3,cex=3) #16)
      title(main=gear.names[g],line=title.line)
   }
   show.block.number(block,dat$t[1],line=2)

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
       F.max=max(Fmort,na.rm=TRUE)
       if (F.max > r)
         F.max = max(F.max,r)
       Fyield = seq(0,F.max,0.01*F.max)
       yield = Fyield*K*(1.0-Fyield/r) # equilibirum yield at F
       print(paste(K,r))
       print(head(cbind(Fyield,yield)))
       obsC = rowSums(dat[,obsC.ndx])
       predC = rowSums(dat[,predC.ndx])
       xrange = c(0,F.max)
       yrange = c(0,max(obsC,predC,yield,na.rm=TRUE))
       Flab = parse(text=paste("Total~Fishing~Mortality~(","y^-1",")",sep=""))
       plot(xrange,yrange,type='n', xlab=Flab, ylab="Total Yield (mt)")

       double.lines(Fmort,predC,bcol="darkgreen",fcol="lightgreen",lwd=5) 
       points(Fmort,obsC,col= "darkgreen",pch=3,cex=2)
       points(Fmort,predC,col="darkgreen",pch=16)
   #   wmaxC = which(obsC==max(obsC,na.rm=TRUE))
   #   lines(c(0,Fmort[wmaxC]),c(0,obsC[wmaxC]),col="red",lty="longdash")
       lines(Fyield,yield,col="red",lwd=3,lty="longdash")

       text(x=Fmort[wy5],y=obsC[wy5],labels=floor(dat$t[wy5]),
             pos=4,offset=0.5,cex=0.8)
       show.block.number(block,0)
   } #if (plot.prod)

   new.devices = devices
   return(new.devices)
}

log.diagnostics=function(file="issams_program.log",ntime=61,dt=1,ngear=5,plot.Fmort=FALSE,plot.prod=FALSE)
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
#  print(length(logsdlogYield))
   sdlogYield = exp(as.numeric(log[logsdlogYield+1]))
#  print(sdlogYield)
   K.pos = grep("^K",log)
   K1 = as.numeric(log[K.pos+2])
   wK1 = which(!is.na(K1))
   K = K1[wK1]
#  print(K)
   r.pos = grep("^r",log)
   r1 = as.numeric(log[r.pos+2])
   wr1 = which(!is.na(r1))
   r = r1[wr1]
#  print(r)
   sdlogQ.pos = grep("^sdlogQ:",log)
   sdlogQ = as.numeric(log[sdlogQ.pos+1])
#  print(paste("sdlogQ",sdlogQ))
  
   max.counter = length(res)
   counter = max.counter
   print(paste(max.counter, "blocks found:"))
#  print(res)

   ncol = (3*ngear+4)
   print(paste("ncol",ncol))
   diag = matrix(nrow=ntime,ncol=ncol)
   print(dim(diag))
   cnames = vector(length=ncol)
   print(cnames)

   c = 'n'
   dev.list = vector(mode="numeric",length=4)
   dev.file.names=c("tmp/est_pop","tmp/est_catch","tmp/est_F")
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
      print(paste("Block:",counter))
   #  print(head(diag))
   #  print(tail(diag))
      print(paste("Displaying block ",counter,sep=""))
      new.devices = plot.diagnostics(as.data.frame(diag),dt=dt,ngear=ngear,
                    sdlogPop=sdlogPop[counter], 
                    sdlogYield=sdlogYield[counter], 
                    sdlogF=sdlogF[counter],
                    sdlogQ=sdlogQ[counter],
                    K=K[counter], r=r[counter],
                    plot.Fmort=plot.Fmort,
                    plot.prod=plot.prod,
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

show.block.number=function(block.number,x,line=3)
{
   mtext(text=paste("(",block.number,")",sep=""),side=1,line=line,
          at=c(x,0),cex=0.8)
}
