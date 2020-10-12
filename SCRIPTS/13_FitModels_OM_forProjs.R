## 13_FitModels_OM_forProjs.R
# HM 200610, helenmoor@gmx.ch

# Fits occupancy models (OMs) to set occupancy at T0/ initialise projections (see SI table S5).
# OMs are fitted: 
# - for mature stands: using the same covariates selected for colonisation probability and Pext ~ 1 
# - for clear-cuts: usign Pcol ~ 1, Pext ~ 1 (only for species tricabie, fomipini, gloesepi, antrseri, phelviti, see SI Fig. S4)



rm(list=ls())

## dirs:------------------------------------------

# for mature stands, load models to select appropriate covariates
INDIR_MOD <- "./MODELS/ForProjections/CEM_mature/"   
INDIR_DAT <- "./DATA/INPUT_DATA/"

OUTDIR <- "./MODELS/ForProjections/OccupancyModels/"

dir.create("MODELS/ForProjections/")
dir.create("MODELS/ForProjections/OccupancyModels/")



### Requirements:------------------------------------------
library(jagsUI) # run jags


### species list:----------------------------------------------------------

modelspecies <- c("amyllapp", "antrseri" ,"fomipini", "fomirose", "gloesepi" , "phelferr", "phelnigr", "phelviti", "phlecent", "tricabie")


### model files (FINAL CE models, to get covariates selected):-------------
myfilenames <- list.files(path = INDIR_MOD, pattern = "_model.RData")

# species
spmodlist <- substr(myfilenames, 0,8)


### LOOP over all models:--------------------------------------------------

for(i in 1:length(myfilenames)){
    
  # current species:
  myspecies= spmodlist[i]

  
  ## load data:----------------------------
  if(myspecies %in% c("phelviti", "gloesepi", "tricabie")){
   minDiam=5
  }else{
    minDiam=10
  }
  load(paste(INDIR_DAT, "ModelInputData_minDiam", minDiam, ".RData" ,sep=""))


  ## Note: for occupancy model, we need an area only offset,
  # contained in offsets as second column ( area scaled to 0.2 ha)
  
  ## prepare Ydata:--------------------------------------------------------

  sel <- which(dimnames(Ydat.mature)[[2]] == myspecies)
  dimnames(Ydat.mature)[[2]] [sel]
  # make a matrix of sites x timepoints for target species: for both mature stands and recently clear-cut stands
  YdatMature <- cbind(Ydat.mature[,sel,1], Ydat.mature[,sel,2])
  YdatClearcuts <- cbind(Ydat.clearcuts[,sel,1], Ydat.clearcuts[,sel,2])
  
  #  Z data, with uncertain 0s (set NA):
  Zdat.in.mature <- YdatMature
  Zdat.in.mature[Zdat.in.mature == 0] <- NA
  # inits for Z, with uncertain 1s (set NA; complementary to Zdat.in):
  Zinits.mature <- YdatMature
  Zinits.mature[Zinits.mature == 1] <- NA
  
  #  Z data, with uncertain 0s (set NA):
  Zdat.in.clearcuts <- YdatClearcuts
  Zdat.in.clearcuts[Zdat.in.clearcuts == 0] <- NA
  # inits for Z, wiht uncertain 1s (set NA; complementary to Zdat.in):
  Zinits.clearcuts <- YdatClearcuts
  Zinits.clearcuts[Zinits.clearcuts == 1] <- NA
  
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

  
  
  #### Set modelling parameters MATURE STANDS ########################################################
  
  ## LOAD FINAL model to get all data needed:--------------------------------------
  load(paste(INDIR_MOD, myfilenames[i], sep=""))  # loads out1
  
  # # Covariates needed for Pcol: (betas)
  colparms <- dimnames(out1$data$XEnvCol)[[2]] #  contains an interecpt
  
  # covariates for occupancy probability:
  XEnvCol <- model.matrix(reformulate(paste("1 + ", colparms[-1], sep="") ), data = Xdat.mature.standardised.Pcol)  
  KEnvCol <- ncol(XEnvCol)
  
  n.cells = length(Xdat.mature.standardised.Pcol[,1])
  
  
  ##### MCMC settings:------------------------------------------
  # ni= timesteps; nt= thinning, nb= burn in; nc = nr of chains
  ni <- 300000   ;   nt <- 100   ;   nb <- 30000  ;  nc <- 3  # final model settings, 3000 samples
  
  # nr of posterior samples:
  n.sims <- ((ni-nb)/nt)*nc
  #n.cells <- length(Ydat.in[,1]) # nr of plots
  
  
  #### Parameters monitored:------------------------------------------
  params <- c( "betan", "n.occ", "detpar") #   "psi"
  ## "n.occ" is the summary metric for replicated simulated data, for PPC 
  


    #### make JAGS DATA
    
    win.data <- list(Y = YdatMature[,2],  ncells = n.cells, offsets = offsets.mature,
                     X = XEnvCol, KEnvCol = KEnvCol,
                     Z = Zdat.in.mature[,2],
                     det=det, ndet=ndet)

    #### INITS 
    inits <- function() list(betan=rep(-0.366, KEnvCol),  detpar=2,
                             Z = Zinits.mature[,2])

    
    ## RUN
    outOM <- jags(win.data, 
             inits= inits, parameters.to.save = params,    
             parallel = TRUE, 
             model.file = paste("SCRIPTS/JAGSscripts/OM_plot_mature_PPC.txt",sep=""),   
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, verbose = T, store.data= TRUE)
    
    # # plot parameters
    # filenamefig <- paste(OUTDIR, "OM_mature_",spmodlist[i], "_parms.jpg", sep="")
    #     jpeg(filenamefig, width = 800, height = 550, units = "px")
    #     myparms <- outOM$parameters[!grepl("n.occ" , outOM$parameters) & !grepl("deviance" , outOM$parameters) ]  # display all but psi values and deviance
    #     whiskerplot(outOM,parameters=myparms)
    #     dev.off()

    # save models            
    save(outOM, file=paste(OUTDIR, "OM_mature_", spmodlist[i], ".RData", sep=""))
    
    # store summary
    sink(file=  paste(OUTDIR,"OM_mature_", myspecies , "_ModelSummary.txt", sep=""))
    print(outOM)
    sink()
    
    
    ### For recently clear-cut stands ############################################################
    # (only for species tricabie, fomipini, gloesepi, antrseri, phelviti)    
        
    if(myspecies %in% c("tricabie", "fomipini", "gloesepi", "antrseri", "phelviti")){    
    # intercept-only models
    n.cells = length(Xdat.clearcuts.standardised.Pcol[,1]) # nr of plots
    
    ##### MCMC settings:------------------------------------------
    # ni= timesteps; nt= thinning, nb= burn in; nc = nr of chains
    ni <- 300000   ;   nt <- 100   ;   nb <- 200000  ;  nc <- 3  # final model settings
    
    # nr of posterior samples:
    n.sims <- ((ni-nb)/nt)*nc


    #### Parameters monitored:------------------------------------------
    params <- c( "betan", "n.occ", "detpar") #   "psi"
    ## "n.occ" is the summary metric for replicated simulated data, for PPC 
   

    #### make JAGS DATA

    win.data <- list(Y = YdatClearcuts[,2],  ncells = n.cells, offsets = offsets.clearcuts,
                     Z = Zdat.in.clearcuts[,2],
                     det=det, ndet=ndet)

    #### INITS 
    inits <- function() list(betan=rep(-0.366, 1),  detpar=2,
                             Z = Zinits.clearcuts[,2])

    
    ## RUN
    outOM <- jags(win.data, 
             inits= inits, parameters.to.save = params,    
             parallel = TRUE, 
             model.file = paste("SCRIPTS/JAGSscripts/OM_plot_clearcuts_PPC.txt",sep=""),   
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, verbose = T, store.data= TRUE)
    
   # # plot parameters 
   #  filenamefig <- paste(OUTDIR, "OM_clearcuts_",spmodlist[i], "_parms.jpg", sep="")
   #      jpeg(filenamefig, width = 800, height = 550, units = "px")
   #      myparms <- outOM$parameters[!grepl("n.occ" , outOM$parameters) & !grepl("deviance" , outOM$parameters) ]  # display all but psi values and deviance
   #      whiskerplot(outOM,parameters=myparms)
   #      dev.off()

    # save model:            
    save(outOM, file=paste(OUTDIR, "OM_clearcuts_", spmodlist[i], ".RData", sep=""))
    
    # store summary:
    sink(file=  paste(OUTDIR,"OM_clearcuts_", myspecies , "_ModelSummary.txt", sep=""))
    print(outOM)
    sink()
    
    } # end clearcut OM fitting
    
    
    
} # end species loop