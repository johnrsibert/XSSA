fit.plot<-function()
{
   graphics.off()
   cnames = c("t","logF","logN1","logN2","logPredC","logObsC")
   initial = read.table("initial.dat")
   colnames(initial) = cnames
   print(head(initial))
   final = read.table("final.dat")
   colnames(final) = cnames
   print(head(final))

   x11()
   plot(final$t,exp(final$logObsC),xlab="Elapsed Quarter",ylab="Catch")
   lines(final$t,exp(final$logPredC),col="blue")

   
   x11()
   plot(initial$t,exp(initial$logF),xlab="Elapsed Quarter",ylab="F",ylim=c(0,1))
   lines(initial$t,exp(final$logF),col="blue")

   x11()
   plot(final$t,(final$logN2),xlab="Elapsed Quarter",ylab="N2",type='l')
   x11()
   plot(final$t,(final$logN1),xlab="Elapsed Quarter",ylab="N1",type='l')
   lines(final$t,(final$logN1),col="blue")

   return(cbind(initial,final))
}
