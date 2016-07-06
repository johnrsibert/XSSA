print(" ",quote=FALSE)
print("##############################################",quote=FALSE)
require(TMB)
#compile("issams.cpp",flags="-O0 -g",safebounds=TRUE,safeunload=TRUE)
compile("issams.cpp",flags="-O3 -g",safebounds=FALSE,safeunload=TRUE)
dyn.load(dynlib("issams"))

#print("data:")
#print(data)
#print("starting parameters:")
#print(paste("number of parameters",length(parameters)))
#print(names(parameters))
#print(parameters)

#source("read_dat_file.R")
issams.dat = read.dat.file()
print(names(issams.dat))
parameters = issams.dat$parameters
data = issams.dat$data

print("MakeADFun-------------------------------------",quote=FALSE)
obj = MakeADFun(data,parameters,random=c("U"),DLL="issams")
#obj = MakeADFun(data,parameters,DLL="issams")

#obj = ADMB.comp(data,"issams")

obj$control=list(trace=1,eval.max=1,iter.max=2)
Objective.function.value = obj$fn(obj$par)

print("starting optimization-------------------------",quote=FALSE)
st=system.time(opt <- nlminb(obj$par,obj$fn,obj$gr)) #,control=obj$control))
#st=system.time(opt <- nlminb(admb.par,obj$fn,obj$gr,control=list(trace=1)))
print("finished optimization-------------------------",quote=FALSE)
print("Timings:")
print(st)
print(paste("Objective function value =",opt$objective))
sd.rep = sdreport(obj)
max.grad = max(sd.rep$gradient.fixed)
std = rbind(summary(sd.rep,select="fixed"),
            summary(sd.rep,select="report"))
print(paste("Number of parameters = ",length(opt$par),
            " Objective function value = ",opt$objective,
            " Maximum gradient component = ",max.grad,sep=""),quote=FALSE)
print(std)
#print(paste("Convergence:",as.logical(opt$convergence)),quote=FALSE)

#print("Starting MCMC --------------------------------",quote=FALSE)
#rwm = mcmc(obj=obj, nsim=50000, algorithm='RWM') 
#, params.init=opt$par, alpha=.08, diagnostic=TRUE)

print.nll.vector=function(obj)
{
   for (i in 1:(obj$report()$nll_count))
   {
      print(paste(signif(obj$report()$nll_vector[i],15),i,sep=" "),quote=FALSE)
   }
}

## 
## 
## Error in if (log(runif(1)) < fn.new - fn.cur) { 
##   missing value where TRUE/FALSE needed
## In addition: There were 50 or more warnings (use warnings() to see the first 50)
## Timing stopped at: 421.661 0.277 421.819 
## > warnings()
## Warning messages:
## 1: In destructive_Chol_update(L, H, t) :
##   Cholmod warning 'matrix not positive definite' at file ../Supernodal/t_cholmod_super_numeric.c, line 729
## 2: In destructive_Chol_update(L, H, t) :
##   Cholmod warning 'matrix not positive definite' at file ../Supernodal/t_cholmod_super_numeric.c, line 729
## 
## 




#save(obj,opt,std,file="issams-fit.Rdata")

#make.diagnostics=function(residuals)
#{
#   resid.names =c("pop","K","forcing","F1","F2","F3","F4","F5","predC1","predC2","predC3","predC4","predC5","obsC1","obsC2","obsC3","obsC4","obsC5")
#   colnames(residuals)=resid.names
#   head(residuals)
#}
