source("/home/jsibert/Projects/xssa/scripts/utils.R")

plot.diagnostics=function(file="diagnostics.dat")
{
   dat = read.table(file=file,header=TRUE)

   if (dat$t[1] == 1)
   {
      dat$t = (1951.875 + dat$t*0.25)
   #  print(head(dat))
   }
   prop = exp(dat$pop1)/(exp(dat$pop1)+exp(dat$pop2))
   print(summary(cbind(dat,prop)))
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
   plot(dat$t,prop,ylim=c(0,1),type='l',col="red",lwd=lwd)

   nice.ts.plot(dat$t,exp(dat$F1),bcol="darkgreen",fcol="lightgreen",lwd=lwd)

   nice.ts.plot(dat$t,cbind(dat$predC1,dat$obsC1),bcol="orange4",fcol="orange",lwd=lwd)


   return(dat)
}
