
#XSSA
##State-space Stock Assessment with "off-line coupling"


Development of state-space stock assessment models for application to open systems where the stock under analysis is connected to a wider population. The preliminary applications is the yellowfin tuna population which supports fisheries in the Main Hawaiian Islands. The MHI population is connected to a larger tuna population for which stock assessments are available. The estimated population size in areas surrounding the MHI used to "force" the population dynamics.

###Models

* issams6 - single population 6 parameter model with index of abundance; MSY, Fmsy parameterization
* xssams - two population model with explicit exchange; __deprecated__


###Directories

* ADMB - source (.tpl files) for ADMB implementation including several __deprecated__ versions
*	TMB - source (.R, .cpp and .dat files) for TMB port of ADMB implmentation
*	scripts - R and bash scripts for plotting and data manipulation
*	Reports - LaTeX and pdf documentaiton of models and results
* dat - data file to run isams6
*	YPR - yield per recruit analysis of Main Hawaiian Islands yellowfing tuna fishery

### To Do
* Modify estimation of Q to allow estimate to exceed 1.0 for augmented catch models
* Rationalize TMB code with issam6.tpl
* Compute Widely Applicable Information Criterion (WAIC)
* Update analysis in the event that WCPFC uptades regional yellowfin assessment


