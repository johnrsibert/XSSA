require(TMB)
compile("xssams.cpp")
dyn.load(dynlib("xssams"))

   get.field = function(p)
   {
     field = sca[p]
     return(field)
   }
   
   get.numeric.field<-function(p)
   {
      ret = as.numeric(get.field(p))
      return(ret)
   }

sca = scan(file="xssams.dat",comment.char="#",what="raw",quiet=TRUE)
catch.file.name = get.field(1)
forcing.file.name = get.field(2)
ntime = get.numeric.field(3)
ngear = get.numeric.field(4)
dt = get.numeric.field(5)

data = structure(list(
ntime = ntime,
ngear = ngear,
dt = dt,
fr = get.numeric.field(6),
#init.log.T12 = get.numeric.field(7),
#phase.logT12 = get.numeric.field(8),
#init.logT21 = get.numeric.field(9),
#phase.logT21 = get.numeric.field(10) ,
#init.logr = get.numeric.field(11),
#phase.logr = get.numeric.field(12),
#init.logK = get.numeric.field(13),
#phase.logK = get.numeric.field(14),
#init.logsdlogF=vector(length=ngear),
#phase.logsdlogF = -1,
#init.logsdlogPop = vector(length=2),
#phase.logsdlogPop = -1,
#init.logsdlogYield=vector(length=ngear),
#phase.logsdlogYield = -1,
#init.LmeanProportion.local = 0.0,
#phase.LmeanProportion.local = -1,
#init.logsdLProportion.local = 0.0,
#phase.logsdLProportion.local = -7,
ObsCatch=matrix(nrow=ngear,ncol=ntime),
ImmigrantBiomass=vector(length=ntime),
maxtime = ntime,
lengthU = ntime*(ngear+2),
ss = 1.0/dt,
# set up U indexing
#Fndxl = vector(length=ntime),
#Fndxu = vector(length=ntime),
Fndxl = seq(0,(ntime-1)*(ngear),ngear)
Fndxu = seq(0,(ntime-1)*(ngear),ngear)+(ngear-1)
utPop1 = ngear*ntime,
utPop2 = ngear*ntime + ntime
))

nn = 14
#for (g in 1:ngear)
#{
   #nn = nn+1
   #data$init.logsdlogF[g] = get.numeric.field(nn)
#}
#nn = nn+1
#data$phase.logsdlogF = get.numeric.field(nn)

#nn = nn+1
#data$init.logsdlogPop[1] = get.numeric.field(nn)
#nn = nn+1
#data$init.logsdlogPop[2] = get.numeric.field(nn)
#nn = nn+1
#data$phase.logsdlogPop = get.numeric.field(nn)
#for (g in 1:ngear)
#{
   #nn = nn+1
   #data$init.logsdlogYield[g] = get.numeric.field(nn)
#}
#nn = nn+1
#data$phase.logsdlogYield = get.numeric.field(nn)
#nn = nn+1
#data$init.LmeanProportionLocal = get.numeric.field(nn)
#nn = nn+1
#data$phase.LmeanProportionLocal = get.numeric.field(nn)
#nn = nn+1
#data$init.logsdLProportion.local = get.numeric.field(nn)
#nn = nn+1
#data$phase.logsdLProportionLocal = get.numeric.field(nn)

for ( t in 1:ntime)
{
   data$Fndxl[t] = (t-1)*ngear+1
   data$Fndxu[t] = t*ngear
}

#print(nn)
#print(data$phase.logsdLProportion.local)

data$ImmigrantBiomass = as.matrix(read.table(file=forcing.file.name))[data$fr,]
data$ObsCatch = as.matrix(read.table(file=catch.file.name))
#rm(list=c("sca","g","ngear","nn","ntime","dt","t"))

parameters = list(
logT12 = -12.0,
logT21 = -12.0,
logr = -2.3026,
logK = 12.207,
logsdlogF = rep(1.1,data$ngear),
logsdlogPop = c(1.2,1.3),
logsdlogYield=rep(1.4,data$ngear),
LmeanProportionLocal = 2.1972,
logsdLProportionLocal = 2.0,
U = vector(length=data$lengthU,mode="numeric")
)
parameters$U=rep(0.0,data$lengthU)

obj <- MakeADFun(data,parameters,random=c("U"),DLL="xssams")
