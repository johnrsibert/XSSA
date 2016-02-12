print(" ",quote=FALSE)
print("##############################################",quote=FALSE)
require(TMB)
compile("issams.cpp",flags="-O0 -g",safebounds=TRUE,safeunload=TRUE)
dyn.load(dynlib("issams"))

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


dat.file = "./issams.dat"
sca = scan(file=dat.file,comment.char="#",what="raw",quiet=TRUE)
print(paste("Read",length(sca),"items from ",dat.file))

data = list()
phases = list()

data$ngear = get.numeric.field()
nobs.gear=data$ngear
data$ntime = get.numeric.field()
ntime=data$ntime
data$dt = get.numeric.field()
print(paste(nobs.gear,ntime,data$dt))
tcatch=matrix(nrow=nobs.gear,ncol=ntime)
print(dim(data$obs.catch))
for (g in 1:nobs.gear)
{
   for (t in 1:ntime)
   {
      tcatch[g,t] = get.numeric.field()
   }
}
data$use.klingons=get.numeric.field()
data$use.klingon.multiplier=get.numeric.field()

nzero = ntime;
ziter = 0;
while (nzero > 0)
{
   ziter = ziter + 1
   nzero = 0;
   for (g in 1:nobs.gear)
      for (t in 2:ntime)
         if ( (tcatch[g,t] <= 0.0) 
              && (tcatch[g,t-1] > 0.0) && (tcatch[g,t+1] > 0.0) )
         {
            nzero = nzero + 1 
            print(paste(nzero,ziter))
            tcatch[g,t] = 0.5*(tcatch[g,t-1] + tcatch[g,t+1])
            print(paste(nzero, " catch for gear ", g ," at time ", t,
                  " set to ", data$obs.catch[g,t],sep=""))
         }
}
print(paste("Zero catch bridging instances:", nzero))
ZeroCatch = 1.0
data$obs_catch = t(log(tcatch+ZeroCatch))
data$first_year = vector(length=nobs.gear)
data$last_year = vector(length=nobs.gear)
for (g in 1:nobs.gear)
{
   data$first_year[g] = 0
   data$last_year[g] = ntime-1
}

forcing.matrix=matrix(nrow=9,ncol=data$ntime)
for (r in 1:9)
{
   for (y in 1:ntime)
   {
      forcing.matrix[r,y] = get.numeric.field()
   }
}
data$fr = get.numeric.field()
data$immigrant_biomass = forcing.matrix[data$fr,]

data$use_mean_forcing = get.numeric.field()
mean.immigrant.biomass = mean(forcing.matrix[data$fr]);
maximum.immigrant.biomass = max(forcing.matrix[data$fr]);
if (data$use_mean_forcing)
   data$immigrant_biomass = mean.immigrant.biomass;
print(data$immigrant.biomass)


data$phase_Fmsy = get.numeric.field()
phases = c(phases,data$phase_Fmsy)
data$init_Fmsy = get.numeric.field()

data$use_r_prior = get.numeric.field()
r.prior = get.numeric.field()
data$logr_prior = log(r.prior)
sdr.prior = get.numeric.field()
data$varr_prior = sdr.prior*sdr.prior

data$phase_MSY = get.numeric.field()
phases = c(phases,data$phase_MSY)
data$init_MSY = get.numeric.field()

data$phase_sdlogProc = get.numeric.field()
phases = c(phases,data$phase_sdlogProc)
data$init_sdlogProc = get.numeric.field()

data$phase_sdlogYield = get.numeric.field()
phases = c(phases,data$phase_sdlogYield)
data$init_sdlogYield = get.numeric.field()

data$use_Q = get.numeric.field()
data$phase_Q = get.numeric.field()
phases = c(phases,data$phase_Q)
data$init_Q = get.numeric.field()

data$use_robustY = get.numeric.field()
phase_pcon = get.numeric.field()
init_pcon = get.numeric.field()
data$pcon = init_pcon

print(paste(field.counter,"input fields processed"))

data$lengthU = ntime*(nobs.gear+1)
# set up U indexing starting at 0
data$Fndxl = seq(0,(ntime-1)*(nobs.gear),nobs.gear)
data$Fndxu = data$Fndxl+(nobs.gear-1) 
data$utPop = data$Fndxu[length(data$Fndxu)]+1

parameters = list(
  logFmsy = log(data$init_Fmsy),
  logMSY = log(data$init_MSY),
  logsdlogProc = log(data$init_sdlogProc),
  logsdlogYield = log(data$init_sdlogYield),
  logQ = log(data$init_Q)
)
parameters$U=rep(0.0,data$lengthU)

r = 2.0*data$init_Fmsy
logK = log(4.0*data$init_MSY/r)
ut = 0;
for (t in 1:ntime)
{   
   for (g in 1:data$ngear)
   {
      ut = ut + 1
      parameters$U[ut] =  -5.0
   }
}
for (t in 1:ntime)
{   
   ut = ut + 1
   parameters$U[ut] = logK
}

print("data:")
print(data)
print("starting parameters:")
print(paste("number of parameters",length(parameters)))
print(names(parameters))
#print(parameters)

print("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",quote=FALSE)
obj = MakeADFun(data,parameters,random=c("U"),DLL="issams")

print("starting optimization^^^^^^^^^^^^^^^^^^^^^^^^^",quote=FALSE)
cont.list=list(trace=1)#,abs.tol=1e-3,rel.tol=1e-3)
opt = nlminb(obj$par,obj$fn,obj$gr,control=cont.list)
print("opt^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",quote=FALSE)
print(opt)


make.diagnostics=function(residuals)
{
   resid.names =c("pop","K","forcing","F1","F2","F3","F4","F5","predC1","predC2","predC3","predC4","predC5","obsC1","obsC2","obsC3","obsC4","obsC5")
   colnames(residuals)=resid.names
   head(residuals)
}
