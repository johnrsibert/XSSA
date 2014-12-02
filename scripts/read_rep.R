read.rep<-function(fit="plot",path=".")
{

   get.fields<-function(sca,nfield=1)
   {
     p<-.field.counter
   # print(paste("p =",p, "; nfield = ",nfield))
     field<-sca[p:(p+nfield-1)]
     p<-p+nfield;
   # print(paste("p =",p))
     .field.counter<<-p
     return(field)
   }
   
   get.numeric.fields<-function(sca,nfield=1)
   {
      ret<-as.numeric(get.fields(sca,nfield))
      return(ret)
   }

  fit.path<-paste(path,.Platform$file.sep,fit,".rep",sep="")
  print(paste("Reading",fit.path),quote=F)
  sca<-scan(file=fit.path,comment.char="#",what="raw",quiet=TRUE)
  print(paste("Scan read",length(sca),"fields"),quote=F)
  .field.counter<<-1
  ret<-list()

  ret$fit<-fit.path

# Number of time periods
  ret$ntime<-get.numeric.fields(sca)

# Year 1
  ret$year1<-get.numeric.fields(sca)

# Number of regions
  ret$nregion<-get.numeric.fields(sca)

# Number of age classes
  ret$nage<-get.numeric.fields(sca)

# Number of recruitments per year
  ret$nrecruit<-get.numeric.fields(sca)

# Number of fisheries
  ret$nfish<-get.numeric.fields(sca)

# Number of realizations per fishery
  ret$fishery.realizations<-vector(length=ret$nfish)
  for (f in 1:ret$nfish)
  {
     ret$fishery.realizations[f]<-get.numeric.fields(sca)
  }

# Region for each fishery
  ret$fishery.regions<-vector(length=ret$nfish)
  for (f in 1:ret$nfish)
  {
     ret$fishery.regions[f]<-get.numeric.fields(sca)
  }

# Time of each realization by fishery (down)
  ret$realization.time<-matrix(nrow=ret$nfish,ncol=ret$ntime)
  for (f in 1:ret$nfish)
  {
     for (ft in 1:ret$fishery.realization[f])
     {
        rt <-get.numeric.fields(sca)
        t <- as.integer(((rt-ret$year1))*4)+1
        ret$realization.time[f,t] <- rt
      # print(paste(f,rt,t,ret$realization.time[f,t]))
     }
  }
 

# Mean lengths at age
  ret$mean.length.at.age<-vector(length=ret$nage)
  for (a in 1:ret$nage)
  {
     ret$mean.length.at.age[a] <- get.numeric.fields(sca)
  }

# SD of length at age
  ret$sd.length.at.age<-vector(length=ret$nage)
  for (a in 1:ret$nage)
  {
     ret$sd.length.at.age[a] <- get.numeric.fields(sca)
  }

# Mean weights at age
  ret$mean.weight.at.age<-vector(length=ret$nage)
  for (a in 1:ret$nage)
  {
     ret$mean.weight.at.age[a] <- get.numeric.fields(sca)
  }

# Natural mortality at age
  ret$M.at.age<-vector(length=ret$nage)
  for (a in 1:ret$nage)
  {
     ret$M.at.age[a] <- get.numeric.fields(sca)
  }

# Selectivity by age class (across) and fishery (down)
  ret$selectivity<-matrix(nrow=ret$nfish,ncol=ret$nage)
  for (f in 1:ret$nfish)
  {
     for (a in 1:ret$nage)
     {
        ret$selectivity[f,a] <- get.numeric.fields(sca)
     }
  } 

# Catchability by realization (across) by fishery (down)
# ??
  ret$catchability<-matrix(nrow=ret$nfish,ncol=ret$ntime)
  for (f in 1:ret$nfish)
  {
     for (ft in 1:ret$fishery.realization[f])
     {
        q <-get.numeric.fields(sca)
        t <- as.integer(((ret$realization.time[f,ft]-ret$year1))*4)+1
        ret$catchability[f,t] <- q
     #  print(paste(f,t,ret$catchability[f,t]))
     }
  }

# Catchability+effort dev. by realization (across) by fishery (down)
# ??
  ret$qE.devs<-matrix(nrow=ret$nfish,ncol=ret$ntime)
  ret$catchability<-matrix(nrow=ret$nfish,ncol=ret$ntime)
  for (f in 1:ret$nfish)
  {
     for (ft in 1:ret$fishery.realization[f])
     {
        qE <-get.numeric.fields(sca)
        t <- as.integer(((ret$realization.time[f,ft]-ret$year1))*4)+1
        ret$qE.devs[f,ft] <- qE
     #  print(paste(f,ft,ret$qE.devs[f,ft]))
     }
   }

# Fishing mortality by age class (across) and year (down)
  ret$F.at.age<-matrix(nrow=ret$ntime,ncol=ret$nage)
  for (y in 1:ret$ntime)
  {
     for (a in 1:ret$nage)
     { 
        fa<-get.numeric.fields(sca)
     #  print(paste(y,a,fa))
        ret$F.at.age[y,a] <-fa
     }
  }  

# Fishing mortality by age class (across), year (down) and region (block)
  nfmort <- ret$nregion*ret$ntime*ret$nage
  F.vect<-as.vector(get.numeric.fields(sca,nfmort))
  ret$F.at.age.and.region <- array(data=F.vect,
                                   dim=c(ret$nage,ret$ntime,ret$nregion),
                                   dimnames=c("age","year","region"))
  # tail(t(Fayr[,,1]))


# Population number by age (across), year (down) and region
  npop <- ret$nregion*ret$ntime*ret$nage
  ret$pop.at.age.and.region<-array(data=get.numeric.fields(sca,npop),
                                   dim=c(ret$nage,ret$ntime,ret$nregion),
                                   dimnames=c("age","year","region"))
  # head(t(junk$pop.at.age.and.region[,,1]))
  # tail(t(junk$pop.at.age.and.region[,,1]))

# Absolute biomass by region (across) and year (down)
# ret$abs.biomass.by.region<-matrix(data=get.numeric.fields(sca,ret$ntime*ret$nregion),
#                                    nrow=ret$ntime,ncol=ret$nregion)

# Exploitable population by fishery (across) and year (down)
  data <- get.numeric.fields(sca,ret$ntime*ret$nfish)
  ret$exploitable.pop<-matrix(data=data,nrow=ret$ntime,ncol=ret$nfish)

# Recruitment
  data <- get.numeric.fields(sca,ret$ntime*ret$nregion)
  ret$recruitment<-matrix(data=data,nrow=ret$ntime,ncol=ret$nregion,byrow=TRUE)

# Total biomass
  data <- get.numeric.fields(sca,ret$ntime*ret$nregion)
  ret$biomass<-matrix(data=data,nrow=ret$ntime,ncol=ret$nregion,byrow=TRUE)

# Adult biomass
  data <- get.numeric.fields(sca,ret$ntime*ret$nregion)
  ret$adult.biomass<-matrix(data=data,nrow=ret$ntime,ncol=ret$nregion,byrow=TRUE)

# Relative biomass by region (across) and year (down)
  data <- get.numeric.fields(sca,ret$ntime*ret$nregion)
  print(tail(data))
  ret$relative.biomass<-matrix(data=data,nrow=ret$ntime,ncol=ret$nregion,byrow=TRUE)

# Observed catch by fishery (down) and time (across)
# ret$obs.catch.by.fishery<-matrix(0,nrow=ret$nfish,ncol=ret$ntime)
# for (f in 1:ret$nfish)
# {
#    for (ft in 1:ret$fishery.realization[f])
#    {
#       dat <-get.numeric.fields(sca)
#       t <- as.integer(((ret$realization.time[f,ft]-ret$year1))*4)+1
#       ret$obs.catch.by.fishery[f,t] <- dat
#    #  print(paste(f,ft,ret$qE.devs[f,ft]))
#    }
#  }
#  tmp<-get.numeric.fields(sca)
#  print(tmp)
  data<-get.numeric.fields(sca,ret$nfish*ret$ntime)
  print(tail(data))
  ret$obs.catch.by.fishery<-matrix(data=data,nrow=ret$nfish,ncol=ret$ntime,byrow=TRUE)
   tmp<-get.numeric.fields(sca)
   print(tmp)

# Predicted catch by fishery (down) and time (across)
# data<-get.numeric.fields(sca,ret$nfish*ret$ntime)
# print(tail(data))
# ret$pred.catch.by.fishery<-matrix(data=data,nrow=ret$nfish,ncol=ret$ntime,byrow=TRUE)


  ############################################################
  print(paste("Final field counter =",.field.counter),quote=F)
  return(ret)
}


# MULTIFAN-CL Viewer
# 2.0
# Frq file = yft.frq
# Input par file = junk
# Output par file = junk

# Observed CPUE by fishery (down) and time (across)
# Predicted CPUE by fishery (down) and time (across)
# Yield analysis option: 0=none, 1=Bev&Holt, 2=Pella Tomlinson
# Beverton-Holt stock-recruitment relationship report
# alpha = 4.533e+07  beta =  2.221e+04  steepness = 9.393e-01
# Observed spawning Biomass
# Observed recruitment
# Spawning Biomass
# Predicted recruitment
# Beverton-Holt yield analysis report
# MSY
# F multiplier at MSY
# F at MSY
# Total biomass at MSY
# Adult biomass at MSY
# Effort multiplier
# Equilibrium yield
# Equilibrium adult biomass
# Equilibrium total biomass
# Adult biomass over adult biomass at MSY
# Total biomass over total biomass at MSY
# Aggregate F over F at MSY
# Aggregate F
# Yield per recruit report
# Effort multiplier
# Yield per recruit
# Tag reporting rates
# Grouping indicator (0 = no grouping, >0 = grouping)
# Time series variation in reporting rates (0 = no, >0 = yes)
# Reporting rates by fishery (no time series variation)
# No. of time periods associated with tag returns
# Time periods associated with grouped tag returns
# Observed tag returns by time period (across) by fishery groupings (down)
# Predicted tag returns by time period (across) by fishery groupings (down)
# Maximum time at liberty
# Observed vs predicted tag returns by time at liberty
# Movement analysis
# Region 1
# Region 2
# Region 3
# Region 4
# Region 5
# Region 6
# Total biomass in absence of fishing
# Adult biomass in absence of fishing
# Exploitable populations in absence of fishing 
# Predicted catch for interaction analysis by fishery (down) and time (across)
