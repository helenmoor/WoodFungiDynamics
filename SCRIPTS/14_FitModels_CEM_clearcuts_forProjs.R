## 14_FitModels_CEM_clearcuts_forProjs.R
# HM 200610, helenmoor@gmx.ch

# Fits colonization-extinction models (CEMs) for recently clear-cut stands used in projections 
# for species tricabie, fomipini, gloesepi, antrseri, phelviti.
# Here, only effects of DW volume and stand age were tested for, see SI table S6.


rm(list=ls())

### Requirements:
library(jagsUI) # run jags

### dirs:--------------------------------------------------------------------

INDIR <- "./DATA/INPUT_DATA/"
OUTDIR <- "./MODELS/ForProjections/CEM_clearcuts/"

dir.create("./MODELS/ForProjections/CEM_clearcuts/")

### get species list:----------------------------------------------------------

modelspecies <- c("amyllapp", "antrseri" ,"fomipini", "fomirose", "gloesepi" , "phelferr", "phelnigr", "phelviti", "phlecent", "tricabie")

# CEMs for recently clear-cut stands could only be fitted for the five most common species:
species2model <- c( "antrseri", "fomipini", "gloesepi", "phelviti","tricabie")
species2modelnum <- which(modelspecies %in% species2model)
modelspecies[species2modelnum]


# covariates----------------------------------------------------------------------------

#colparmslist - (still listed  for all species, but only fitting models for 5 species)

colparmslist <- list( c(),                               # amylapp
                      c(),                               # antrseri
                      c("DWVol_m3perHa_T2_SpruceConif"), # fomipini
                      c(),                               # fomirose
                      c("DWVol_m3perHa_T2_SpruceConif"), # gloesepi
                      c(),                               # phelferr
                      c(),                               # phelnigr
                      c("AgeT2"),                        # phelviti
                      c(),                               # phlecent
                      c())                               # tricabie



extparmslist <- list( c(),                               # amylapp
                      c(),                               # antrseri
                      c(),                               # fomipini
                      c(),                               # fomirose
                      c(),                               # gloesepi
                      c(),                               # phelferr
                      c(),                               # phelnigr
                      c(),                               # phelviti
                      c(),                               # phlecent
                      c())                               # tricabie



### loop over species:--------------------------------------------------------------
for(myloop in species2modelnum){ 
    
    myspecies= modelspecies[myloop]
    
    # for phelviti, "gloesepi", "tricabie" use minDiam = 5 cm dataset, all other species minDiam = 10cm
    if(myspecies %in% c("phelviti", "gloesepi", "tricabie")){
        minDiam = 5
    }else{
        minDiam = 10}
    # load data:------------------------------------
    load(file=paste(INDIR, "ModelInputData_minDiam", minDiam, ".RData", sep=""))
    
    
    ## prepare Ydata:--------------------------------------------------------
    sel <- which(dimnames(Ydat.clearcuts)[[2]] == myspecies)
    dimnames(Ydat.clearcuts)[[2]] [sel]
    # make a matrix of sites x timepoints for target species:
    Ydat <- cbind(Ydat.clearcuts[,sel,1], Ydat.clearcuts[,sel,2])
    
    ## Model older stands only:-----------------------------------
    # covariates:
    # rename for indataset
    XdatscaledPcol.in <- Xdat.clearcuts.standardised.Pcol
    XdatscaledPext.in <- Xdat.clearcuts.standardised.Pext  
    # indata for Y: observed
    Ydat.in <- Ydat   
    # indata for Z: uncertain zeros (set NA)
    Zdat.in <- Ydat.in                                              
    Zdat.in[Zdat.in == 0] <- NA
    # Z inits, uncertain 1s (set NA ; complementary to Zdat.in):
    Zinits <- Ydat.in
    Zinits[Zinits == 1] <- NA

    # offsets
    offsets.in <- offsets.clearcuts
 
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
  
  
    ##### MCMC settings:------------------------------------------
    # ni= timesteps; nt= thinning, nb= burn in; nc = nr of chains
    ni <- 300000   ;   nt <- 10   ;   nb <- 200000  ;  nc <- 3  # final model settings, 15000 samples
    
    # nr of posterior samples:
    n.sims <- ((ni-nb)/nt)*nc
    n.cells <- length(Ydat.in[,1]) # nr of plots
    
    
    #### Parameters monitored:------------------------------------------
    params <- c( "beta", "gamma", "lomega", "n.occT1", "n.occT2", "n.ext", "n.col","detpar") 
    ## "n.occT1", "n.occT2", "n.ext", "n.col" are summary metrix for replicated simulated data, for PPC
    
    ### set required covariates - pick from lists:----------------------

    colparms <-  colparmslist[[myloop]] 
    extparms <-  extparmslist[[myloop]] 
    
    
    ## run Jags:---------------------------
    
    ## Prepare indata:
    # Colonization:
    XEnvCol <- model.matrix(reformulate(paste(c("1", colparms), sep="") ), data = XdatscaledPcol.in)   
    KEnvCol <- ncol(XEnvCol)
    # Extinction
    XEnvExt <- model.matrix(reformulate(paste(c("1",  extparms), sep="") ), data = XdatscaledPext.in) 
    KEnvExt <- ncol(XEnvExt)
    
    #### make JAGS DATA
    win.data <- list(Y = Ydat.in,  ncells = dim(Ydat.in)[1], offsets = offsets.in,
                     XEnvCol = XEnvCol, KEnvCol = KEnvCol, XEnvExt = XEnvExt, KEnvExt = KEnvExt,
                     Z = Zdat.in,
                     det=det, ndet=ndet)
    
    #### INITS
    inits <- function() list(beta=rep(-0.366, KEnvCol), gamma=rep(-2, KEnvExt), omega = 0, detpar=2,
                             Z = Zinits)
    
    
    ## different JAGS scripts depending on model type:
    
    ## Pext ~ 1
    ## Intercept only on Pext
    if(length(extparms)==0){
        
        ## Pext ~ 1, Pcol ~ 1:
        if(length(colparms)==0){
            
            ##### RUN
            out1 <- jags(win.data, 
                         inits= inits, parameters.to.save = params,    
                         parallel = TRUE, 
                         model.file = paste("SCRIPTS/JAGSscripts/CEdet_Plot_ColExtIntercept_PPC.txt",sep=""),   
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, verbose = T,
                         store.data=TRUE)
            
        }else if(length(colparms)>0){
            ## Pext ~ 1, Pcol ~ covs: 
            
            ##### RUN
            out1 <- jags(win.data, 
                         inits= inits, parameters.to.save = params,    
                         parallel = TRUE, 
                         model.file = paste("SCRIPTS/JAGSscripts/CEdet_Plot_ExtIntercept_PPC.txt",sep=""),   
                         n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, verbose = T,
                         store.data=TRUE)
        } # end if(length(colparms)>0)
        
    }# end model fitting with Pext ~ 1
 
    
    # # plot parameter estimates:---------------------------
    # filenamefig <- paste(OUTDIR,myspecies , "_parms.jpg", sep="")
    # jpeg(filenamefig, width = 800, height = 550, units = "px")
    # myparms <- out1$parameters[!grepl("psi", out1$parameters) & !grepl("deviance", out1$parameters)  & !grepl("n.occT1" , out1$parameters)  & !grepl("n.occT2" , out1$parameters)  & !grepl("n.ext" , out1$parameters)  & !grepl("n.col", out1$parameters) ]  # display all but psi values and deviance
    # whiskerplot(out1,parameters=myparms)
    # dev.off()
 
    ## save model
    save(out1, file=  paste(OUTDIR,myspecies , "_model.RData", sep="") )
    
    # store summary:
    sink(file=  paste(OUTDIR,myspecies , "_ModelSummary.txt", sep=""))
    print(out1)
    sink()
    
    rm(out1)
    
} # end species loop

