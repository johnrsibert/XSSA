require(TMB)
compile("xssams.cpp")
dyn.load(dynlib("xssams"))

field.counter <<- 0

   get.field = function()
   {
     field.counter <<- field.counter + 1
     field = sca[field.counter]
   # print(paste(field.counter,field))
     return(field)
   }
   
   get.numeric.field<-function()
   {
      ret = as.numeric(get.field())
      return(ret)
   }

logit<-function(p)
{
   return(log(p/(1-p)))
}

alogit<-function(alpha)
{
   return(1/(1+exp(-alpha)))
}


dat.file = "../run/xssams.dat"
sca = scan(file=dat.file,comment.char="#",what="raw",quiet=TRUE)
print(paste("Read",length(sca),"items from ",dat.file))
data = list()

data$ngear = get.numeric.field()
ngear=data$ngear
data$ntime = get.numeric.field()
ntime=data$ntime
data$dt = get.numeric.field()
data$obs.catch=matrix(nrow=ngear,ncol=ntime)
for (g in 1:ngear)
{
   for (y in 1:ntime)
   {
      data$obs.catch[g,y] = get.numeric.field()
   }
}
nzero = ntime;
ziter = 0;
while (nzero > 0)
{
   ziter = ziter + 1
   nzero = 0;
   for (g in 1:ngear)
      for (t in 2:ntime)
         if ( (data$obs.catch[g,t] <= 0.0) 
              && (data$obs.catch[g,t-1] > 0.0) && (data$obs.catch[g,t+1] > 0.0) )
         {
            nzero = nzero + 1 
            print(paste(nzero,ziter))
            data$obs.catch[g,t] = 0.5*(data$obs.catch[g,t-1] + data$obs.catch[f,t+1])
            print(paste(nzero, " catch for gear ", g ," at time ", t,
                  " set to ", data$obs.catch[g,t],sep=""))
         }
}
print(paste("Zero catch bridging instances:", nzero))
ZeroCatch = 1.0
data$obs.catch = log(data$obs.catch+ZeroCatch);


forcing.matrix=matrix(nrow=9,ncol=data$ntime)
for (r in 1:9)
{
   for (y in 1:ntime)
   {
      forcing.matrix[r,y] = get.numeric.field()
   }
}
data$fr = get.numeric.field()
data$immigrant.biomass = forcing.matrix[data$fr,]

data$use.mean.forcing = get.numeric.field()
mean.immigrant.biomass = mean(forcing.matrix[data$fr]);
maximum.immigrant.biomass = max(forcing.matrix[data$fr]);
if (data$use.mean.forcing)
   data$immigrant.biomass = mean.immigrant.biomass;
#print(data$immigrant.biomass)

data$phase.T12 = get.numeric.field()
data$init.T12 = get.numeric.field()
data$phase.T21 = get.numeric.field() 
data$init.T21 = get.numeric.field()
data$phase.r = get.numeric.field()
data$init.r = get.numeric.field()
data$phase.K = get.numeric.field()
data$init.K = get.numeric.field()


data$phase.sdlogF = get.numeric.field()
data$init.sdlogF = get.numeric.field()
data$phase.sdlogPop = get.numeric.field()
data$init.sdlogPop = get.numeric.field()
data$phase.sdlogYield = get.numeric.field()
data$init.sdlogYield = get.numeric.field()
data$phase.meanProportion.local = get.numeric.field()
data$init.meanProportion.local =  get.numeric.field()
data$phase.sdProportion.local = get.numeric.field()
data$init.sdProportion.local = get.numeric.field()
data$phase.qProp = get.numeric.field()
data$init.qProp = get.numeric.field()
data$use.robustY = get.numeric.field()
data$phase.pfat = get.numeric.field()
data$init.pfat = vector(length=ngear)
for (g in 1:ngear)
{
#  print(g)
   data$init.pfat[g] = get.numeric.field()
}
print(paste(field.counter,"input fields processed"))

data$maxtime = ntime
data$lengthU = ntime*(ngear+2)
data$ss = 1.0/data$dt
# set up U indexing starting at 0
data$Fndxl = seq(0,(ntime-1)*(ngear),ngear)
data$Fndxu = data$Fndxl+(ngear-1) 
data$utPop1 = ngear*ntime - 1
data$utPop2 = data$utPop1 + ntime


parameters = list()
parameters$logT12 = log(data$init.T12+1e-10)
parameters$logT21 = log(data$init.T21+1e-10)
parameters$logr = log(data$init.r)
parameters$logK = log(data$init.K)
parameters$logsdlogF = log(data$init.sdlogF)
parameters$logsdlogYield = log(data$init.sdlogYield)
parameters$logsdlogPop = log(data$init.sdlogPop)
parameters$LmeanProportionLocal = logit(data$init.meanProportion.local)
parameters$logsdLProportionLocal = log(logit(data$init.sdLProportion.local))
parameters$qProp = data$init.qProp
if (!data$use.robustY)
{
   data$phase.pfat = -1;
   data$init.pfat = 1e-25
}
parameters$Lpfat = logit(data$init.pfat)

parameters$U = vector(length=data$lengthU,mode="numeric")

parameters$U=rep(0.0,data$lengthU)

obj = MakeADFun(data,parameters,random=c("U"),DLL="xssams")

#obj <- MakeADFun(data,parameters,random=c("U"),DLL="sam")
print(names(obj))
print(obj$par)
lower <- obj$par*0-Inf
upper <- obj$par*0+Inf
print(paste(lower,upper))
#lower["rho"] <- 0.01
#upper["rho"] <- 0.99

system.time(opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper))
#system.time(opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper))
rep<-sdreport(obj)
rep
