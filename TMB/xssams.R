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

dat.file = "xssams.dat"
sca = scan(file=dat.file,comment.char="#",what="raw",quiet=TRUE)
print(paste("Read",length(sca),"items from ",dat.file))
data = list()

catch.file.name = get.field()
forcing.file.name = get.field()

data$ntime = get.numeric.field()
data$ngear = get.numeric.field()
data$dt = get.numeric.field()
data$fr = get.numeric.field()
data$init.logT12 = get.numeric.field()
data$phase.logT12 = get.numeric.field()
data$init.logT21 = get.numeric.field()
data$phase.logT21 = get.numeric.field() 
data$init.logr = get.numeric.field()
data$phase.logr = get.numeric.field()
data$init.logK = get.numeric.field()
data$phase.logK = get.numeric.field()

ntime=data$ntime
ngear=data$ngear

data$init.logsdlogF=vector(length=ngear)
for (g in 1:ngear)
   data$init.logsdlogF[g] = get.numeric.field()
data$phase.logsdlogF = get.numeric.field()
data$init.logsdlogPop = get.numeric.field()
data$phase.logsdlogPop = get.numeric.field()
data$init.logsdlogYield=vector(length=ngear)
for (g in 1:ngear)
   data$init.logsdlogYield[g] = get.numeric.field()
data$phase.logsdlogYield = get.numeric.field()
data$init.LmeanProportion.local =  get.numeric.field()
data$phase.LmeanProportion.local = get.numeric.field()
data$init.logsdLProportion.local = get.numeric.field()
data$phase.logsdLProportion.local = get.numeric.field()
data$use.robustF = get.numeric.field()
data$init.Lpfat = vector(length=ngear)
for (g in 1:ngear)
   data$init.Lpfat[g] = get.numeric.field()
data$phase.pfat = get.numeric.field()
print(paste(field.counter,"input fields processed"))

data$ImmigrantBiomass = as.matrix(read.table(file=forcing.file.name))[data$fr,]
data$ObsCatch = as.matrix(read.table(file=catch.file.name))

data$maxtime = ntime
data$lengthU = ntime*(ngear+2)
data$ss = 1.0/data$dt
# set up U indexing starting at 0
data$Fndxl = seq(0,(ntime-1)*(ngear),ngear)
data$Fndxu = data$Fndxl+(ngear-1) 
data$utPop1 = ngear*ntime - 1
data$utPop2 = data$utPop1 + ntime


parameters = list()
parameters$logT12 = data$init.logT12
parameters$logT21 = data$init.logT21
parameters$logr = data$init.logr
parameters$logK = data$init.logK
parameters$logsdlogF = data$init.logsdlogF
parameters$logsdlogPop = data$init.logsdlogPop
parameters$logsdlogYield = data$init.logsdlogYield
parameters$LmeanProportionLocal = data$init.LmeanProportion.local
parameters$logsdLProportionLocal = data$init.logsdLProportion.local
parameters$Lpfat = data$init.Lpfat
parameters$U = vector(length=data$lengthU,mode="numeric")

parameters$U=rep(0.0,data$lengthU)

obj = MakeADFun(data,parameters,random=c("U"),DLL="xssams")
