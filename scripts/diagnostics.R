source("/home/jsibert/Projects/xssa/scripts/utils.R")

plot.diagnostics=function(dat=NULL,file="diagnostics.dat")
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
   old.par = par(no.readonly = TRUE) 
   par(mar=c(3,4.5,0,0)+0.1)
   np = 3
   lm = layout(matrix(c(1:np),ncol=1,byrow=TRUE))
   layout.show(lm)
   lwd = 3

   nice.ts.plot(dat$t,cbind(dat$pop1,dat$pop2),bcol="blue",fcol="lightblue",lwd=lwd)
   par("new"=TRUE)
   plot(dat$t,dat$propL,ylim=c(0,1),type='l',col="red",lwd=lwd)

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
            diag[t,c] = as.numeric(log[fc])
         }
      }
      print(paste("Block:",counter))
      print(head(diag))
      print(tail(diag))
      plot.diagnostics(as.data.frame(diag))
      title(main=paste(log[afters[counter]+1],"entries"))
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



