source("/home/jsibert/Projects/xssa/scripts/utils.R")

plot.diagnostics=function(dat=NULL,file="diagnostics.dat",ngear)
{
   if (is.null(dat))
     dat = read.table(file=file,header=TRUE)

   if (dat$t[1] == 1)
   {
      dat$t = (1951.875 + dat$t*0.25)
      print(names(dat))
   }
#  dat$pop1 = exp(dat$pop1)
#  dat$pop2 = exp(dat$pop2)

   pop.dev = 0
   catch.dev = 0
   F.dev = 0
   lwd = 3

   width = 9.0
   height = 11.0
   if (pop.dev > 0)
      dev.set(pop.dev)
   x11(width=width,height=height)
   pop.dev = dev.cur()
   print(paste("Current device",dev.cur()))
   print(paste("pop device",pop.dev))
   old.par = par(no.readonly = TRUE) 
   par(mar=c(3,4.5,0,0)+0.1)
   np = 3
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)

   nice.ts.plot(dat$t,dat$pop1,lwd=lwd)
   title(main="Population 1")
   nice.ts.plot(dat$t,dat$pop2,lwd=lwd)
   title(main="Population 2")
   nice.ts.plot(dat$t,(dat$pop1+dat$pop2),lwd=lwd)
   par("new"=TRUE)
   plot(dat$t,dat$propL,type='l',ann=FALSE,axes=FALSE,
        col="red",lwd=2,lty="dotted")
   title(main="Total Population")
   save.png.plot("est_pop",width=width,height=height)

   width = 9.0
   height =11.0
   if (catch.dev > 0)
      dev.set(catch.dev)
   x11(width=width,height=height)
   catch.dev = dev.cur()

   par(mar=c(3,4.5,0,0)+0.1)
   np = ngear
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   for (g in 1:ngear)
   {
      nice.ts.plot(dat$t,dat[,(g+9)],bcol="darkgreen",fcol="lightgreen",lwd=lwd)
      points(dat$t,dat[,(g+14)],col= "darkgreen",pch=16)
      title(main=paste("Catch, gear",g))
   }
   save.png.plot("est_catch",width=width,height=height)

   width = 9.0
   height =11.0
   if (F.dev > 0)
      dev.set(F.dev)
   x11(width=width,height=height)
   F.dev = dev.cur()
   par(mar=c(3,4.5,0,0)+0.1)
   np = ngear
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   for (g in 1:ngear)
   {
      nice.ts.plot(dat$t,dat[,(g+4)],bcol="orange4",fcol="orange",lwd=lwd)
      title(main=paste("F mort, gear",g))
   }
   save.png.plot("est_F",width=width,height=height)
#  return(dat)
#  return(cbind(dat$pop1,dat$pop2))
}

plot.diagnostics.1=function(dat=NULL,file="diagnostics.dat")
{
   if (is.null(dat))
     dat = read.table(file=file,header=TRUE)

   if (dat$t[1] == 1)
   {
      dat$t = (1951.875 + dat$t*0.25)
   #  print(head(dat))
   }
#  prop = exp(dat$pop1)/(exp(dat$pop1)+exp(dat$pop2))
#  print(summary(cbind(dat,prop)))
   dat$pop1 = exp(dat$pop1)
#  w = which(dat$pop1 <= 2)
#  print(w)
#  is.na(dat$pop1[w])
   dat$pop2 = exp(dat$pop2)
#  print(head(cbind(dat$pop1,dat$pop2))) 
#  dat$predC1=exp(dat$predC1)
#  dat$obsC1=exp(dat$obsC1)

   width = 6.5
   height = 9.0
   x11(width=width,height=height)
   print(paste("Current device",dev.cur()))
   old.par = par(no.readonly = TRUE) 
   par(mar=c(3,4.5,0,0)+0.1)
   np = 3
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   lwd = 3

   nice.ts.plot(dat$t,cbind(dat$pop1,dat$pop2),bcol="blue",fcol="lightblue",lwd=lwd)
   par("new"=TRUE)
   plot(dat$t,dat$propL,ylim=c(0,1),type='l',col="red",lwd=lwd,axes=FALSE,ann=FALSE)

   nice.ts.plot(dat$t,exp(dat$F1),bcol="darkgreen",fcol="lightgreen",lwd=lwd)

   nice.ts.plot(dat$t,cbind(dat$predC1,dat$obsC1),bcol="orange4",fcol="orange",lwd=lwd)

   return(dat)
#  return(cbind(dat$pop1,dat$pop2))
}

log.diagnostics=function(file="xssams_program.log",ntime=244,ngear=5)
{
      
   print(paste("Scanning file",file))
   log = scan(file,what="character")
   res = grep("Residuals:",log)
   afters = grep("after",log)

   counter = 1
   max.counter = length(res)
   print(paste(max.counter, "blocks found:"))
   print(res)

   ncol = (3*ngear+4)
   diag = matrix(nrow=ntime,ncol=ncol)
   cnames = vector(length=ncol)


   c = 'n'
   while (c != 'q')
   {
      fc = res[counter]
      for (c in 1:ncol)
      {
         fc = fc + 1
      #  print(paste(fc,log[fc])) 
         cnames[c] = log[fc]
      }
      print(cnames)
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
      print(head(diag))
      print(tail(diag))
      plot.diagnostics(as.data.frame(diag),ngear=ngear)
#     title(main=paste(log[afters[counter]+1],"entries"))
      c = readline("Next, back or quit (enter n,b, or q):")
      print(paste(c," entered"))
      if (c == 'n')
         counter = counter + 1
      else if (c == 'b')
         counter = counter - 1
   
      if (counter < 1)
         counter = max.counter
      else if (counter > max.counter)
         counter = 1
   } #while
   graphics.off()
#  return(counter)
}



