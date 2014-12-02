plot.sel<-function()
{
   nage<-28
   nfleet<-17
   sel<-read.table("sel.dat",head=F)
   head(sel)
   plot(c(1,28),c(0,1),type='n',xlab="Age (quarter)",ylab="Selectivity")
   for (i in 1:nfleet)
   {
     lines(seq(1,28,1),sel[i,])

   }
   return(sel)
}
