## 11_FitModels_CEM_main.R
# HM 200610, helenmoor@gmx.ch

# Fits the main colonization-extinction models (CEMs) for all species, including effects of local connectivity for amyllapp and phlecent. See SI table S7 for parameter estimates.

rm(list=ls())

### Requirements:
library(jagsUI) # run jags


### dirs:--------------------------------------------------------------------

INDIR <- "./DATA/INPUT_DATA/"
OUTDIR <- "./MODELS/CEM_main/"

dir.create("MODELS/")
dir.create("MODELS/CEM_main/")


### species list:----------------------------------------------------------

modelspecies <- c("amyllapp", "antrseri" ,"fomipini", "fomirose", "gloesepi" , "phelferr",
"phelnigr", "phelviti", "phlecent", "tricabie")

## LocalConnectivity (here called presence in "cell6") was important only for two species: "amyllapp", "phlecent" (added as separate dataset: Xdat.mature.LocalConnectivity)


#colparmslist: list of covariates used in each model for colonisation probability, by species
colparmslist <- list( c("DWVol_m3perHa_T2_SpruceConif"), # amylapp, plus local connectivity
                   c("DWVol_m3perHa_T2_SpruceConif"), # antrseri
                   c("DWVol_m3perHa_T2_SpruceConif"), # fomipini
                   c("DWVol_m3perHa_T2_SpruceConif", "AgeT2"), # fomirose
                   c("SprucevolMax100m_2013_m3ha"), # gloesepi
                   c("DWVol_m3perHa_T2_SpruceConif", "AgeT2"), # phelferr
                   c("DWVol_m3perHa_T2_SpruceConif", "AgeT2"), # phelnigr
                   c("DWVol_m3perHa_T2_SpruceConif", "AgeT2"), # phelviti
                   c("DWVol_m3perHa_T2_SpruceConif"), # phlecent, plus local connectivity
                   c("DWVol_m3perHa_T2_SpruceConif")) # tricabie

# extparmslist: list of covariates used in each model for extinction probability, by species
extparmslist <- list( c(), # amylapp
                   c(), # antrseri
                   c(), # fomipini
                   c(), # fomirose
                   c(), # gloesepi
                   c(), # phelferr
                   c("DWVol_m3perHa_T2_SpruceConif"), # phelnigr
                   c("DWVol_m3perHa_T2_SpruceConif"), # phelviti
                   c(), # phlecent
                   c("DWVol_m3perHa_T2_SpruceConif","MeanDecaySwe_T2")) # tricabie


### loop over species:--------------------------------------------------------------

for(myloop in 1:length(modelspecies)){
  
  myspecies <- modelspecies[myloop] # current species being modelled
  
  # for gloesepi, phelviti, tricabie use minDiam = 5 cm dataset, all other species minDiam = 10cm
  if(myspecies %in% c("phelviti", "tricabie", "gloesepi")){
   minDiam = 5
    }else{
    minDiam = 10}

  # load data:------------------------------------
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
    # indata for X: rename Xdat for indataset
    XdatscaledPcol.in <- Xdat.mature.standardised.Pcol
    XdatscaledPext.in <- Xdat.mature.standardised.Pext 
    
    # indata for Y (observed occurrences): observed
    Ydat.in <- Ydat     
    
    # indata for Z (true unobserved occurrences), where 0s are uncertain (set NA):
    Zdat.in <- Ydat
    Zdat.in[Zdat.in == 0] <- NA

    # offsets:
    offsets.in <- offsets.mature
    
    ## prepare LocalConnectivity (Cell6) data and inits:------------------------------------
    PresCell6Dat <- as.data.frame(Xdat.mature.LocalConnectivity[,myspecies])
    

    
    ##### MCMC settings:------------------------------------------
    # ni= timesteps; nt= thinning, nb= burn in; nc = nr of chains
    ni <- 300000   ;   nt <- 100   ;   nb <- 200000  ;  nc <- 3  # final model settings

    # nr of posterior samples:
    n.sims <- ((ni-nb)/nt)*nc
    n.cells <- length(Ydat.in[,1]) # nr of plots
    

    
    ### set required covariates - pick from lists:----------------------
    # excluding LocalConnectivity (added separately for amyllapp and phlecent)
    colparms <-  colparmslist[[myloop]] # 
    extparms <-  extparmslist[[myloop]] #  
    

        
    #### Parameters monitored:------------------------------------------
    ## "n.occT1", "n.occT2", "n.ext", "n.col" are summary metrix for replicated simulated data, for PPC
    
   if(myspecies %in% c("amyllapp", "phlecent")){
    params <- c( "beta","betaCell6", "gamma", "delta", "lomega", "n.occT1", "n.occT2", "n.ext", "n.col", "detpar") # only "amyllapp", "phlecent" require betaCell6 and delta parameter
   }else{ 
    params <- c( "beta","gamma", "lomega", "n.occT1", "n.occT2", "n.ext", "n.col", "detpar")}
    
    
    
    ## run Jags:---------------------------
    
    ## with covariates on Pext (requires different JAGS script):-----------
    if(length(extparms)>0){
    # (not the case for amyllapp and phlecent)
    # fits models for tricabie, phelviti, phelnigr
        
      # Fit model:
      XEnvCol <- model.matrix(reformulate(paste(c("1", colparms), sep="") ), data = XdatscaledPcol.in)   
      KEnvCol <- ncol(XEnvCol)
      # # Extinction (default: intercept only)
      XEnvExt <- model.matrix(reformulate(paste("1 +",  extparms, sep="") ), data = XdatscaledPext.in) 
      KEnvExt <- ncol(XEnvExt)
      
    #### make JAGS DATA
    win.data <- list(Y = Ydat.in,  ncells = dim(Ydat.in)[1], offsets = offsets.in,
                       XEnvCol = XEnvCol, KEnvCol = KEnvCol, XEnvExt = XEnvExt, KEnvExt = KEnvExt,
                       Z = Zdat.in,
                       det=det, ndet=ndet)
 
    #### INITS: initial values for parameters 
    # inits for Z, where 1s are uncertain (set NA; complementary to Zdat.in):
    Zinits <- Ydat
    Zinits[Zinits == 1] <- NA
    # inits for parameters
    inits <- function() list(beta=rep(-0.366, KEnvCol), gamma=rep(-2, KEnvExt), omega = 0, detpar=2,
                               delta=0,
                               Z = Zinits)
      
    #### RUN
      
      out1 <- jags(win.data, 
                   inits= inits, parameters.to.save = params,    
                   parallel = TRUE, 
                   model.file = paste("SCRIPTS/JAGSscripts/CEdet_Plot_PPC.txt",sep=""),   
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, verbose = T,
                   store.data=TRUE)
      
      # finish model fitting with Pext ~ covs
      
    }else if(length(extparms)==0){  
      ##  in case Pext is modelled as intercept only (requires different JAGS script):-----------
      
        # fit models for amyllapp and phlecent (add LocalConnectivity)
        if(myspecies %in% c("amyllapp", "phlecent")){
            
         
        # Fit model
        # Colonisation:
        XEnvCol <- model.matrix(reformulate(paste(c("1", colparms), sep="") ), data = XdatscaledPcol.in)   
        KEnvCol <- ncol(XEnvCol)
        # Extinction (LocalConnectivity (PresCell6) is added as separate dataset)
        XEnvExt <- model.matrix(~ 1 , data = XdatscaledPext.in)   
        KEnvExt <- ncol(XEnvExt)
      
        #### make JAGS DATA
        win.data <- list(Y = Ydat.in,  ncells = dim(Ydat.in)[1], offsets = offsets.in,
                       XEnvCol = XEnvCol, KEnvCol = KEnvCol, XEnvExt = XEnvExt, KEnvExt = KEnvExt,
                       PresCell6 = PresCell6Dat[,1],
                       Z = Zdat.in,
                       det=det, ndet=ndet)
 
        #### INITS 
        # inits for Z, where 1s are uncertain (set NA; complementary to Zdat.in):
        Zinits <- Ydat
        Zinits[Zinits == 1] <- NA
        # inits for LocalConn (Presence in Cell6): complementary to data (set NA where there is data; set 0 where data is NA)
        PresCell6inits <-  rep(NA, times=length(PresCell6Dat[,1]))
        PresCell6inits[which(is.na(PresCell6Dat[,1]))] <-  0
        
        inits <- function() list(beta=rep(-0.366, KEnvCol), gamma=rep(-2, KEnvExt), omega = 0, detpar=2,
                               betaCell6 = rnorm(0,1),
                               delta=0,
                               Z = Zinits,
                               PresCell6 = PresCell6inits)     

      
        #### RUN (CEdet_Plot_ExtIntercept_PPC_wLocalConn.txt)
      
        out1 <- jags(win.data, 
                   inits= inits, parameters.to.save = params,    
                   parallel = TRUE, 
                   model.file = paste("SCRIPTS/JAGSscripts/CEdet_Plot_ExtIntercept_PPC_wLocalConn.txt",sep=""),   
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, verbose = T,
                   store.data=TRUE)
        }else{
        
        # Fit models with Pext ~ 1 and no effect of LocalConnectivity: 
        # for fomipini, gloesepi, antrseri, phelferr, fomirose
            
        # Fit model
        # Colonisation:
        XEnvCol <- model.matrix(reformulate(paste(c("1", colparms), sep="") ), data = XdatscaledPcol.in)   
        KEnvCol <- ncol(XEnvCol)
        # Extinction: (intercept only)
        XEnvExt <- model.matrix(~ 1 , data = XdatscaledPext.in)   
        KEnvExt <- ncol(XEnvExt)
      
        #### make JAGS DATA
        win.data <- list(Y = Ydat.in,  ncells = dim(Ydat.in)[1], offsets = offsets.in,
                       XEnvCol = XEnvCol, KEnvCol = KEnvCol, XEnvExt = XEnvExt, KEnvExt = KEnvExt,
                       Z = Zdat.in,
                       det=det, ndet=ndet)
 
        #### INITS 
        # inits for Z, where 1s are uncertain (set NA; complementary to Zdat.in):
        Zinits <- Ydat
        Zinits[Zinits == 1] <- NA
        inits <- function() list(beta=rep(-0.366, KEnvCol), gamma=rep(-2, KEnvExt), omega = 0, detpar=2,
                               delta=0,
                               Z = Zinits)     

      
        #### RUN (CEdet_Plot_ExtIntercept_PPC.txt)
      
        out1 <- jags(win.data, 
                   inits= inits, parameters.to.save = params,    
                   parallel = TRUE, 
                   model.file = paste("SCRIPTS/JAGSscripts/CEdet_Plot_ExtIntercept_PPC.txt",sep=""),   
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, verbose = T,
                   store.data=TRUE)    
            
        } # end fit models with Pext ~ 1 and no effect of LocalConnectivity: 
        
        
    } # end model fitting with Pext ~ 1
    

    # # plot parameter estimates:---------------------------
    # filenamefig <- paste(OUTDIR,myspecies , "_parms.jpg", sep="")
    # jpeg(filenamefig, width = 800, height = 550, units = "px")
    # myparms <- out1$parameters[!grepl("psi", out1$parameters) & !grepl("deviance", out1$parameters)  & !grepl("n.occT1" , out1$parameters)  & !grepl("n.occT2" , out1$parameters)  & !grepl("n.ext" , out1$parameters)  & !grepl("n.col", out1$parameters) ]  # display all but psi values, deviance, and summary metrics
    # whiskerplot(out1,parameters=myparms)
    # dev.off()

    
    ## save model
    save(out1, file=  paste(OUTDIR,myspecies , "_model.RData", sep="") )
      
    # save summary: 
    sink(file=  paste(OUTDIR, myspecies , "_ModelSummary.txt", sep=""))
    print(out1)
    sink()
    
    rm(out1)
  
} # end species loop
  
                      