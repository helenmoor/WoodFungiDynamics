## 10_FitModels_BaseRates.R
# HM 200610, helenmoor@gmx.ch

# This script loops over all 10 model species, loads the required data (minimum diameter depends on species), and fits an intercept-only model (base colonization and extinction rates). See Table 1 for parameter estimates.

rm(list=ls())

### Requirements:
library(jagsUI) # run jags


### dirs:--------------------------------------------------------------------

INDIR <- "./DATA/INPUT_DATA/"
OUTDIR <- "./MODELS/BaseRates/"

dir.create("MODELS/")
dir.create("MODELS/BaseRates/")

### species list:----------------------------------------------------------

modelspecies <- c("amyllapp", "antrseri" ,"fomipini", "fomirose", "gloesepi" , "phelferr",
"phelnigr", "phelviti", "phlecent", "tricabie")


### Model setup:-----------------------------------------------------------
# Manually set required final model covariates for each species.
# Base rate models are intercept only models; colonization parameters and extinction parameters remain empty.

#colparmslist: list of covariates used in each model for colonisation probability, by species
colparmslist <- list( c(), # amylapp
                   c(), # antrseri
                   c(), # fomipini
                   c(), # fomirose
                   c(), # gloesepi
                   c(), # phelferr
                   c(), # phelnigr
                   c(), # phelviti
                   c(), # phlecent
                   c()) # tricabie

# extparmslist: list of covariates used in each model for extinction probability, by species
extparmslist <- list( c(), # amylapp
                   c(), # antrseri
                   c(), # fomipini
                   c(), # fomirose
                   c(), # gloesepi
                   c(), # phelferr
                   c(), # phelnigr
                   c(), # phelviti
                   c(), # phlecent
                   c()) # tricabie




### loop over species:--------------------------------------------------------------
for(myloop in 1:length(modelspecies)){ 
  
  myspecies <- modelspecies[myloop] # current species being modelled
  
  # for gloesepi, phelviti, tricabie use minDiam = 5 cm dataset, all other species minDiam = 10cm
  if(myspecies %in% c("phelviti", "tricabie", "gloesepi")){
   minDiam = 5
    }else{
    minDiam = 10}
  
  ## load data:------------------------------------
  # see DataOverview.html for description of dataset
  load(file=paste(INDIR, "ModelInputData_minDiam", minDiam, ".RData", sep=""))
  
  
  ## detection data:--------------------
    # use Variable "Detected":
      det <- detdat$Detected
      ndet = length(det)
  
    # for tricabie, fomipini and antrseri: hand it species-specific detection data:
      if(myspecies == "tricabie"){
        det <- det[detdat$Species=="tricabie"]
         ndet <- length(det)
      }
      # do the same for fomipini and antrseri:
      if(myspecies == "fomipini"){
        det <- det[detdat$Species=="fomipini"]
         ndet <- length(det)
      }
      if(myspecies == "antrseri"){
        det <- det[detdat$Species=="antrseri"]
         ndet <- length(det)
      }
      # for all other species, keep pooled detection data.
  
  
  ## prepare Ydata:--------------------------------------------------------
  sel <- which(dimnames(Ydat.mature)[[2]] == myspecies)
  dimnames(Ydat.mature)[[2]] [sel]
  # make a matrix of sites x timepoints for target species:
  Ydat <- cbind(Ydat.mature[,sel,1], Ydat.mature[,sel,2])

  
  ## Indata to model:-----------------------------------
    # rename Xdat for indataset
    XdatscaledPcol.in <- Xdat.mature.standardised.Pcol
    XdatscaledPext.in <- Xdat.mature.standardised.Pext 
    
    # indata for Y (observed occurrences): observed
    Ydat.in <- Ydat     
    
    # indata for Z (true unobserved occurrences), where 0s are uncertain (set NA):
    Zdat.in <- Ydat
    Zdat.in[Zdat.in == 0] <- NA

    # offsets:
    offsets.in <- offsets.mature
    

    ##### MCMC settings:------------------------------------------
    # ni= timesteps (incl burnin); nt= thinning, nb= burn in; nc = nr of chains
    ni <- 300000   ;   nt <- 100   ;   nb <- 200000  ;  nc <- 3  # final model settings

    # nr of posterior samples:
    n.sims <- ((ni-nb)/nt)*nc
    n.cells <- length(Ydat.in[,1]) # nr of plots
    
    
    #### Parameters monitored:------------------------------------------
    params <- c( "beta", "gamma", "lomega", "n.occT1", "n.occT2", "n.ext", "n.col", "detpar") 
    # not required: "psi",  "pcol", "pext"
    # "n.occT1", "n.occT2", "n.ext", "n.col" are summary metrics for replicated simulated data, for PPC
    
    ### set required covariates for current species - pick from lists:----------------------

    #names(Xdat.mature.standardised.Pcol)
    colparms <-  colparmslist[[myloop]] 
    extparms <-  extparmslist[[myloop]] 
    
    
    ## run Jags:---------------------------
  
      # Model to fit:
      # Colonisation:
      XEnvCol <- model.matrix(~ 1 , data = XdatscaledPcol.in)   
      KEnvCol <- ncol(XEnvCol) # K parameters to be estimated, including intercept
      # # Extinction:
      XEnvExt <- model.matrix(~ 1 , data = XdatscaledPext.in)   
      KEnvExt <- ncol(XEnvExt)  # K parameters to be estimated, including intercept
      
      #### make JAGS DATA:

      win.data <- list(Y = Ydat.in,  ncells = dim(Ydat.in)[1], offsets = offsets.in,
                       XEnvCol = XEnvCol, KEnvCol = KEnvCol, XEnvExt = XEnvExt, KEnvExt = KEnvExt,
                       Z = Zdat.in,
                       det=det, ndet=ndet)
      
      #### INITS: initial values for parameters 
      
      # inits for Z, where 1s are uncertain (set NA; complementary to Zdat.in):
      Zinits <- Ydat
      Zinits[Zinits == 1] <- NA
      # inits for parameters
      inits <- function() list(beta=-0.366, gamma=-2, omega = 0, detpar=2,
                               Z = Zinits)
      

      
      ##### RUN:
      
      out1 <- jags(win.data, 
                   inits= inits, parameters.to.save = params,    
                   parallel = TRUE, 
                   model.file = paste("SCRIPTS/JAGSscripts/CEdet_Plot_ColExtIntercept_PPC.txt",sep=""),   
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, verbose = T,
                   store.data=TRUE)
    
    # result:
    out1$summary
        
    ## save model:
    save(out1, file=  paste(OUTDIR,myspecies , "_model.RData", sep="") )
      
    # save summary: 
    sink(file=  paste(OUTDIR, myspecies , "_ModelSummary.txt", sep=""))
    print(out1)
    sink()
      
    
    rm(out1)
  
} # end species loop
  