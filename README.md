
#XSSA
##State-space Stock Assessment with population exchange


Development of state-space stock assessment models for application to open systems where the stock under analysis is connected to a wider population. The preliminary applications is the yellowfin tuna population which supports fisheries in the Main Hawaiian Islands. The MHI population is connected to a larger tuna population for which stock assessments are available. The estimated population size in areas surrounding the MHI used to "force" the population dynamics.

###Models

* issams - single population model with index of abundance; MSY, Fmsy parameterization
* issams-dev - single population model with index of abundance; alternative parameteization (K, d_1)
* xssams - two population model with explicit exchange; __deprecated__


###Directories

* ADMB - source (.tpl files) for ADMB implementation
*	TMB - source (.R, .cpp and .dat files) for TMB port of ADMB implmentation
*	scripts - R and bash scripts for plotting and data manipulation
*	Reports - LaTeX and pdf documentaiton of models and results
* run-issams 	TMB code runs to convergence 	21 days ago
*	YPR - yield per recruite of Main Hawaiian Islands yellowfing tuna fishery
*	run-issams - run directory (.dat file) for issams ADMB implementation
*	run-issams-dev - run directory (.dat file) for issams-dev ADMB implementation
*	run-xssams - run directory (.dat file) for xssams ADMB implementation
*	run-xssams-dev - run directory (.dat file) for xssams-dev ADMB implementation`
*	run-issams-msy - __deprecated__
*	run-issams-proc - __deprecated__

