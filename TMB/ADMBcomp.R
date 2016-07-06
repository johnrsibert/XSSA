print(" ",quote=FALSE)
print("################# ADMBcomp.R ####################",quote=FALSE)

require(TMB)
#compile("issams.cpp",flags="-O0 -g",safebounds=TRUE,safeunload=TRUE)
compile("issams.cpp",flags="-O3 -g",safebounds=FALSE,safeunload=TRUE)
dyn.load(dynlib("issams"))

load(file="ADMBfit.Rdata")
#print(names(admb.fit))
admb.parameters = list()
admb.parameters$logFmsy = admb.fit$est[which(admb.fit$names == "logFmsy")]
admb.parameters$logMSY = admb.fit$est[which(admb.fit$names == "logMSY")]
admb.parameters$logsdlogProc = admb.fit$est[which(admb.fit$names == "logsdlogProc")]
admb.parameters$logsdlogYield = admb.fit$est[which(admb.fit$names == "logsdlogYield")]
admb.parameters$logQ = admb.fit$est[which(admb.fit$names == "logQ")]
admb.parameters$U = admb.fit$est[which(admb.fit$names == "U")]
print(names(admb.parameters))

issams.dat = read.dat.file()
#print(names(issams.dat))
#print(names(issams.dat$parameters))

parameters = admb.parameters
data = issams.dat$data

print("MakeADFun-------------------------------------",quote=FALSE)
obj = MakeADFun(data,parameters,random=c("U"),DLL="issams")
#obj = MakeADFun(data,parameters,DLL="issams")

#obj$control=list(trace=1,eval.max=1,iter.max=2)
print("Compute one function evaluation---------------",quote=FALSE)
Objective.function.value = obj$fn(obj$par)
print(paste("Objective function value =",Objective.function.value),quote=FALSE )

print("starting optimization-------------------------",quote=FALSE)
st=system.time(opt <- nlminb(obj$par,obj$fn,obj$gr)) #,control=obj$control))
#st=system.time(opt <- nlminb(admb.par,obj$fn,obj$gr,control=list(trace=1)))
print("finished optimization-------------------------",quote=FALSE)
print("Timings:")
print(st)
print(paste("Objective function value =",Objective.function.value, ", prior to optimization"),quote=FALSE )
print(paste("Objective function value =",opt$objective),quote=FALSE)
sd.rep = sdreport(obj)
max.grad = max(sd.rep$gradient.fixed)
std = rbind(summary(sd.rep,select="fixed"),
            summary(sd.rep,select="report"))
print(paste("ADMB: Number of parameters = ",admb.fit$nopar,
            " Objective function value = ",admb.fit$nlogl,
            " Maximum gradient component = ",admb.fit$maxgrad,sep=""),quote=FALSE)
print(paste(" TMB: Number of parameters = ",length(opt$par),
            " Objective function value = ",opt$objective,
            " Maximum gradient component = ",max.grad,sep=""),quote=FALSE)
print(std)
#print(paste("Convergence:",as.logical(opt$convergence)),quote=FALSE)


