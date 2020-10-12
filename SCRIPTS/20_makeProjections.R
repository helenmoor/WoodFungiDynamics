## 20_makeProjections.R
# HM 200610, helenmoor@gmx.ch

## Make projections of species col-ext dynamics and resulting occupancy for scenario data ExampleStands.
# (results shown in Fig. 5)

## Young production stands (Young Spruce Stands, YSS, 21-63 yrs):
# - common species with occurrences in YSS: antrseri fomipini  gloesepi  tricabie: are modelled with the mature stand model (OM and ColExt) (for rest assume fixed probs)
# - all other species: Pocc = 0 in YSS

## Recently clear-cut stands:
# - separate clearcut stand OM and ColExt models for: antrseri, fomipini, gloesepi, phelviti, tricabie

rm(list=ls())

## packages:---------------------------------------------------------
library(jagsUI) # may not be needed if loading posteriors directly


## dirs:--------------------------------------------------------------

# directories of final fitted models: (model files copied over from PolyporeCEM)
CEM_mature_dir <- "./MODELS/ForProjections/CEM_mature/"
CEM_clearcuts_dir <- "./MODELS/ForProjections/CEM_clearcuts/"
OM_dir <-  "./MODELS/ForProjections/OccupancyModels/"

# projection covariate data dir
PROJDAT_DIR <- "./DATA/SCENARIO_DATA/"

# outdir for results:
OUTdir <- "./PROJECTIONS/"
dir.create("./PROJECTIONS")


# fcts:-------------------------------------------------------------
logit <- function(x) log(x/(1-x))
ilogit <- function(x) exp(x)/(1+exp(x))
cloglog <- function(x) log(-log(1-x))
icloglog <- function(x) 1-exp(-exp(x))


## loop over required dataset

for(minDiam in c(5,10)){  

### INDATA models:----------------------------------------------------------

myspecies <- c("amyllapp", "antrseri", "fomipini", "fomirose", "phelferr", "phelnigr", "phlecent")
  
if(minDiam==5){
  myspecies <- c("gloesepi", "phelviti", "tricabie")
  
}

  
# CEM models (with posteriors) mature stands:
mycemfiles <- list.files(CEM_mature_dir, pattern="_model")
# subset to only target species:
mycemfiles <- mycemfiles[which(grepl(paste(myspecies, collapse="|"), mycemfiles))]

## add clearcuts models: (6 species)
mycemfiles.clearcuts <- list.files(CEM_clearcuts_dir, pattern="_model")
# subset myfiles to only target species:
mycemfiles.clearcuts <- mycemfiles.clearcuts[which(grepl(paste(myspecies, collapse="|"), mycemfiles.clearcuts))]



# OM posteriors:
myomfiles <- list.files(OM_dir, pattern="RData")
# subset myfiles to only target species:
myomfiles <- myomfiles[which(grepl(paste(myspecies, collapse="|"), myomfiles))]
# split into clearcuts and mature OMs:
myomfiles_mature <- myomfiles[grep("mature", myomfiles)]
myomfiles_clearcuts <- myomfiles[grep("clearcuts", myomfiles)]
rm(myomfiles)


### INDATA scenarios:----------------------------------------------------------
# single dataset: 

myscens = "ExampleStands"
myscenfiles <- list.files(paste(PROJDAT_DIR, sep=""), pattern="RData")
nscens <- length(myscens) # 1

  
### INCORPORATE MECHANISTIC ASSUMPTIONS:--------------------

### Mechanistic assumption 0: where deadwoodvol=0 (woodpres=0), Pocc = 0

### Mechanistic assumption 1: age gap treatment depends on species (see above)


#### Projections:-----------------------------------------------------------------

# Load scenario covariate data:--------------------------------------------------

# single scenario (ExampleStands):
scen=1
print(myscens[scen])

## load scenario data: standdat and covdat
load(paste(PROJDAT_DIR,"/", myscenfiles[scen], sep=""))

# covdat contains DWVol for both min diameters, 
# remove the DWVol column that is not right for he current minDiam:
if(minDiam == 10){
    covdat <- covdat[,-3,]
}else if(minDiam == 5){
    covdat <- covdat[,-2,]
}
dimnames(covdat)[[2]][2] <- "DWVol_m3perHa_T2_SpruceConif"  # generalise column name again

# standdat: df with info on stands
# covdat: array with [nstands, covs: interecpt (1), DwVol, StandAge), 11 periods]
# covs are named in same way as in fitted models

### Transform and standardise covariate data:---------------------
## Keep raw covdat to decide on:
# - is there DW or not?
# - ID of mature, recently clear-cut and young production stands in each period

## Standardisation:
## transform, center and scale covs for projections, based on FIN covariate means and sds (stored in in COVinfo), load means and sds:
load(paste("./DATA/INPUT_DATA/COVinfo_minDiam", minDiam, ".RData", sep=""))

## For spruce vol (but not DWVol), add a tiny amount to all occasions where they are zero (to enable log transform):
sel<- which(covdat[,4,]==0, arr.ind = T) 
covdat[sel[,1], 4, sel[,2]] <- 0.001
  
#  make cov datasets: log-transform, center and scale (separate for Pcol and Pext, mature and clearcuts),
## this will results in 4 diff covariate datasets:
# for Pcol and Pext (different log transforms) x for mature and clearcuts (diff means and sds)

## Covs for Pcol and OM Pocc: log-transform all:------------------------

# Covariates for Pcol, mature stands:
covdat.col.mature <- covdat
covdat.col.mature[,2:4,] <- log(covdat.col.mature[,2:4,]) # log-transform all covs 
## center and scale using the COV means and sds from the Finnish data:
sel <- match(dimnames(covdat)[[2]], names(COVinfo$COVinfo.Xdat.mature.Pcol$COVmeans))  # use match for arranging cols
for(i in 2:4){
    COVmean <- COVinfo$COVinfo.Xdat.mature.Pcol$COVmeans[sel[i]]
    COVsd <- COVinfo$COVinfo.Xdat.mature.Pcol$COVsds[sel[i]]
    covdat.col.mature[,i,] <-  (covdat.col.mature[,i,] - COVmean) / COVsd
}
## NOTE: this results in -Inf for DW where DW = 0
# -> use the raw data covdat to decide in each period which stands have no dead wood and set Pocc = 0 there

# Covariates for Pcol, clearcut stands:
covdat.col.clearcuts <- covdat # needed for occupancy in clearcuts stands and clearcuts stand models
covdat.col.clearcuts[,2:4,] <- log(covdat.col.clearcuts[,2:4,]) # log-transform all covs 
# center and scale using the COV means and sds from the Finnish data:
sel <- match(dimnames(covdat)[[2]], names(COVinfo$COVinfo.Xdat.clearcuts.Pcol$COVmeans))  # use match for arranging cols
for(i in 2:4){
    COVmean <- COVinfo$COVinfo.Xdat.clearcuts.Pcol$COVmeans[sel[i]]
    COVsd <- COVinfo$COVinfo.Xdat.clearcuts.Pcol$COVsds[sel[i]]
    covdat.col.clearcuts[,i,] <-  (covdat.col.clearcuts[,i,] - COVmean) / COVsd
}

## Covariates for Pext: log-transform only DWVol (and conn):-----------------

# Covariates for Pext, mature stands:
covdat.ext.mature <- covdat
covdat.ext.mature[,2,] <- log(covdat.ext.mature[,2,]) 
# center and scale using the COV means and sds from the Finnish data:
sel <- match(dimnames(covdat)[[2]], names(COVinfo$COVinfo.Xdat.mature.Pext$COVmeans))  # use match for arranging cols
for(i in 2:3){
    COVmean <- COVinfo$COVinfo.Xdat.mature.Pext$COVmeans[sel[i]]
    COVsd <- COVinfo$COVinfo.Xdat.mature.Pext$COVsds[sel[i]]
    covdat.ext.mature[,i,] <-  (covdat.ext.mature[,i,] - COVmean) / COVsd
}

# Covariates for Pext, clearcuts stands: 
covdat.ext.clearcuts <- covdat
covdat.ext.clearcuts[,2,] <- log(covdat.ext.clearcuts[,2,])
# center and scale using the COV means and sds from the Finnish data:
sel <- match(dimnames(covdat)[[2]], names(COVinfo$COVinfo.Xdat.clearcuts.Pext$COVmeans))  # use match for arranging cols
for(i in 2:3){
    COVmean <- COVinfo$COVinfo.Xdat.clearcuts.Pext$COVmeans[sel[i]]
    COVsd <- COVinfo$COVinfo.Xdat.clearcuts.Pext$COVsds[sel[i]]
    covdat.ext.clearcuts[,i,] <-  (covdat.ext.clearcuts[,i,] - COVmean) / COVsd
}


nstands <- length(standdat$Description)
nperiods <- dim(covdat)[3]



# species loop:--------------------------------------------
for(sp in 1:length(myspecies)){
    
    # doing species:
    species = myspecies[sp]
    print(species)
    
    ## load models and extract posteriors:-----------------------------------------
    # to facilitate matrix multiplication, add columns of zeros for covs not needed
    
    # CEMs:--------------------
    # load model: (object out1)
    load(paste(CEM_mature_dir, mycemfiles[grep(myspecies[sp], mycemfiles)], sep="")) 
    
    # make pcol.posteriors, pext.posteriors (column names are covs, incl intercept)
    # get posteriors from all three chains:
    post <- as.data.frame(out1$samples[[1]])
    post <- rbind(post, out1$samples[[2]])
    post <- rbind(post, out1$samples[[3]])
    
    # if posteriors contain >1000 samples, thin here to 1000:
    if(dim(post)[[1]] >1000){
        set.seed(39) # reproducible
        mysel <- sample(3000,size = 1000, replace = F) # sample joint posterior distribution 
        post <- post[mysel,]
        rm(mysel)
    }
    
    # Covariates needed for Pcol: (betas)
    colparms <- dimnames(out1$data$XEnvCol)[[2]]
    colparms
    # NOTE: if intercept only, colparms is empty - add name for Intercept:
    if(length(colparms) == 0){
        colparms <- "(Intercept)"
    }
    ncolparms <- length(colparms) # number of colparms, including intercept
    
    # Covariates needed for Pext: (gammas)
    extparms <- dimnames(out1$data$XEnvExt)[[2]]
    extparms
    # NOTE: if intercept only, extparms is empty - add name for Intercept:
    if(length(extparms) == 0){
        extparms <- "(Intercept)"
    }
    nextparms <- length(extparms) # including intercept
    
    ## split post into 2 dataframes, and name columns by covariate :
    pcol.posteriors <- as.data.frame(post[,1:ncolparms])#keep as df in case of single column
    names(pcol.posteriors) <- colparms
    
    pext.posteriors <- as.data.frame(post[,(ncolparms+1):(ncolparms+nextparms)]) #keep as df in case of single column
    names(pext.posteriors) <- extparms
    
    ## covariates needed:
    colcovs <- colparms
    extcovs <- extparms
    
    # clean up and make some space:
    rm(post, colparms, extparms)
    rm(out1); gc()
    
    # add zero columns for parameters for covs that are not needed:
    if(length(colcovs) < length(dimnames(covdat)[[2]])){
      missing <- which(!(dimnames(covdat)[[2]] %in% colcovs))
      mydf <- data.frame(matrix(0, ncol = length(missing), nrow = dim(pcol.posteriors)[1]))
      colnames(mydf) <- dimnames(covdat)[[2]][missing]
      # add to pcol.posteriors:
      pcol.posteriors <- cbind(pcol.posteriors, mydf)
      rm(mydf)
      ## reorder to fixed order:
      pcol.posteriors.mature <- pcol.posteriors[,c("(Intercept)" ,"DWVol_m3perHa_T2_SpruceConif",  "AgeT2" ,"SprucevolMax100m_2013_m3ha")]
    }
    # add zero cols for parameters for covs that are not needed:
    if(length(extcovs) < length(dimnames(covdat)[[2]])){
      missing <- which(!(dimnames(covdat)[[2]] %in% extcovs))
      mydf <- data.frame(matrix(0, ncol = length(missing), nrow = dim(pext.posteriors)[1]))
      colnames(mydf) <- dimnames(covdat)[[2]][missing]
      # add to pext.posteriors:
      pext.posteriors <- cbind(pext.posteriors, mydf)
      rm(mydf)
      ## reorder to fixed order:
       pext.posteriors.mature <- pext.posteriors[,c("(Intercept)" ,"DWVol_m3perHa_T2_SpruceConif",  "AgeT2","SprucevolMax100m_2013_m3ha" )]
    }
    
    rm(pcol.posteriors, pext.posteriors)
    
    ## clearcuts model posteriors:----------------
    # -> pcol.posteriors, pext.posteriors (column names are covs, incl intercept)

    if(myspecies[sp] %in% c("antrseri", "fomipini", "gloesepi", "phelviti", "tricabie")){ 
    
        load(paste(CEM_clearcuts_dir, mycemfiles.clearcuts[grep(myspecies[sp] ,mycemfiles.clearcuts)], sep=""))
        
        # make pcol.posteriors, pext.posteriors (column names are covs, incl intercept)
        # get posteriors from all three chains:
        post <- as.data.frame(out1$samples[[1]])
        post <- rbind(post, out1$samples[[2]])
        post <- rbind(post, out1$samples[[3]])
        
        # if posteriors contain >1000 samples, thin here to 1000:
        if(dim(post)[[1]] >1000){
            set.seed(39) # reproducible
            mysel <- sample(3000,size = 1000, replace = F) # sample joint posterior distribution 
            post <- post[mysel,]
            rm(mysel)
        }
        
        # Covariates needed for Pcol: (betas)
        colparms <- dimnames(out1$data$XEnvCol)[[2]]
        colparms
        # NOTE: if intercept only, colparms is empty - add name for Intercept:
        if(length(colparms) == 0){
            colparms <- "(Intercept)"
        }
        ncolparms <- length(colparms)#including intercept
        
        # Covariates needed for Pext: (gammas)
        extparms <- dimnames(out1$data$XEnvExt)[[2]]
        extparms
        # NOTE: if intercept only, extparms is empty - add name for Intercept:
        if(length(extparms) == 0){
            extparms <- "(Intercept)"
        }
        nextparms <- length(extparms) # including intercept
        
        ## split post into 2 dataframes, and name columns by covariate :
        pcol.posteriors <- as.data.frame(post[,1:ncolparms])#keep as df in case single column
        names(pcol.posteriors) <- colparms
        pext.posteriors <- as.data.frame(post[,(ncolparms+1):(ncolparms+nextparms)]) #keep as df in case single column
        names(pext.posteriors) <- extparms
        
        # covariates needed for clearcut CEMs:    
        colcovsclearcuts <- colparms
        extcovsclearcuts <- extparms
    
        rm(out1, post ,colparms, extparms); gc()
    
        # add zero cols for parameters for covs that are not needed:
        if(length(colcovsclearcuts) < length(dimnames(covdat)[[2]])){
        missing <- which(!(dimnames(covdat)[[2]] %in% colcovsclearcuts))
        mydf <- data.frame(matrix(0, ncol = length(missing), nrow = dim(pcol.posteriors)[1]))
        colnames(mydf) <- dimnames(covdat)[[2]][missing]
        # add to pext.posteriors:
        pcol.posteriors <- cbind(pcol.posteriors, mydf)
        rm(mydf)
        ## reorder: 
        pcol.posteriors.clearcuts <- pcol.posteriors[,c("(Intercept)" ,"DWVol_m3perHa_T2_SpruceConif",  "AgeT2" ,"SprucevolMax100m_2013_m3ha")]
        }
        # add zero cols for parameters for covs that are not needed:
        if(length(extcovsclearcuts) < length(dimnames(covdat)[[2]])){
        missing <- which(!(dimnames(covdat)[[2]] %in% extcovsclearcuts))
        mydf <- data.frame(matrix(0, ncol = length(missing), nrow = dim(pext.posteriors)[1]))
        colnames(mydf) <- dimnames(covdat)[[2]][missing]
        # add to pext.posteriors:
        pext.posteriors <- cbind(pext.posteriors, mydf)
        rm(mydf)
        ##  reorder:
        pext.posteriors.clearcuts <- pext.posteriors[,c("(Intercept)" ,"DWVol_m3perHa_T2_SpruceConif",  "AgeT2","SprucevolMax100m_2013_m3ha" )]
        }
        
        rm(pcol.posteriors, pext.posteriors); gc()
        
    } # end load clearcuts model where available
    

    
    ## Occupancy models:-------------------------------------------------------------------
    # OM mature: -> pocc.posteriors (same covariates as Pcol)
    load(paste(OM_dir, myomfiles_mature[grep(myspecies[sp], myomfiles_mature)], sep=""))
    # get posteriors from all three chains:
    post <- as.data.frame(outOM$samples[[1]])
    post <- rbind(post, outOM$samples[[2]])
    post <- rbind(post, outOM$samples[[3]])
    
    # if posteriors contain >1000 samples, thin here to 1000:
    if(dim(post)[[1]] >1000){
        set.seed(13) # reproducible
        mysel <- sample(3000,size = 1000, replace = F) # sample joint posterior distribution 
        post <- post[mysel,]
        rm(mysel)
    }
    
    # Covariates needed for Pocc: (betas)
    occcovs <- dimnames(outOM$data$X)[[2]]
    occcovs
#    identical(occcovs, colcovs)
    # NOTE: if intercept only, add name for Intercept:
    if(length(occcovs) == 0){
        occcovs <- "(Intercept)"
    }
    noccparms <- length(occcovs) # including intercept
    
    # name columns by covariate :
    pocc.posteriors.mature <- as.data.frame(post[,1:noccparms])#keep as df in case single column
    names(pocc.posteriors.mature) <- occcovs
    rm(post)
    
     # IF not all covs needed, add zero cols: (same as colcovs)
    if(length(occcovs) < length(dimnames(covdat)[[2]])){
      missing <- which(!(dimnames(covdat)[[2]] %in% colcovs))
      mydf <- data.frame(matrix(0, ncol = length(missing), nrow = dim(pocc.posteriors.mature)[1]))
      colnames(mydf) <- dimnames(covdat)[[2]][missing]
      # add to pext.posteriors:
      pocc.posteriors.mature <- cbind(pocc.posteriors.mature, mydf)
      rm(mydf)
      ## reorder
       pocc.posteriors.mature <- pocc.posteriors.mature[,c("(Intercept)" ,"DWVol_m3perHa_T2_SpruceConif",  "AgeT2","SprucevolMax100m_2013_m3ha" )]
    }
       
    rm(outOM); gc()
    
    ## OM clearcuts : -> Intercept only
    if(myspecies[sp] %in% c("antrseri", "fomipini", "gloesepi", "phelviti", "tricabie")){ 
    load(paste(OM_dir, myomfiles_clearcuts[grep(myspecies[sp], myomfiles_clearcuts)], sep=""))
    
    # get posteriors from all three chains:
    post <- as.data.frame(outOM$samples[[1]])
    post <- rbind(post, outOM$samples[[2]])
    post <- rbind(post, outOM$samples[[3]])
    
    # if posteriors contain >1000 samples, thin here to 1000:
    if(dim(post)[[1]] >1000){
        set.seed(15) # reproducible
        mysel <- sample(3000,size = 1000, replace = F) # sample joint posterior distribution 
        post <- post[mysel,]
        rm(mysel)
    }
    
    # Intercept only for clearcut OMs:
    # name columns by covariate :
    pocc.posteriors.clearcuts <- as.data.frame(post[,1])#keep as df in case single column
    names(pocc.posteriors.clearcuts) <- "(Intercept)"
    rm(post)

    # add zero cols: (same as colcovs)
    missing <- which(!(dimnames(covdat)[[2]] %in% names(pocc.posteriors.clearcuts)))
    mydf <- data.frame(matrix(0, ncol = length(missing), nrow = dim(pocc.posteriors.clearcuts)[1]))
    colnames(mydf) <- dimnames(covdat)[[2]][missing]
    # add to pext.posteriors:
    pocc.posteriors.clearcuts <- cbind(pocc.posteriors.clearcuts, mydf)
    rm(mydf)

    rm(outOM); gc()
    }
    
    ## all posteriors now contain columns for all covariates in covdat (zeros where cov not needed), 
    # -> this enables to simply use matrix multiplication in each case
    
    
    ## set up results array:-----------------
    nsamp <- dim(pcol.posteriors.mature)[1]
    # array for resulting Probability of Occurrence
    ProbOcc <- array(0,dim=c(nstands,nperiods,nsamp))  #  
    # collect also Pcol
    ProbCol <- array(0,dim=c(nstands,nperiods,nsamp))  #
    # collect also Pext
    ProbExt <- array(1,dim=c(nstands,nperiods,nsamp))   # set default to 1
    
  # identical(dimnames(covdat)[[1]], standdat$Description) # same order
    dimnames(ProbOcc)[[1]] <- standdat$Description  # add stand names to result
    dimnames(ProbCol)[[1]] <- standdat$Description # 
    dimnames(ProbExt)[[1]] <- standdat$Description

    
    ### SIMULATE:--------------------------------------------------
    
    ## PERIOD 1: set initial state with occupancy models:------------------
    # Mature stands: >= 64 years : pocc.posteriors.mature x covs
    # Clearcut stands: 0-20 years: pocc.posteriors.clearcuts x covs (Intercept only)
    # Age gap: >20-63 years: Pocc = 0, except for antrseri, fomipini, gloesepi, tricabie: here use mature stand OM

    ## Clearcut stands: formula: icloglog(intercept)
    # Pocc remains zero, except for c("antrseri", "fomipini", "gloesepi","phelviti", "tricabie")
    if(myspecies[sp] %in% c("antrseri", "fomipini", "gloesepi", "phelviti", "tricabie")){ # then use OM_clearcuts
        # select based on raw covdat:
        sel <- as.numeric(which(covdat[,"AgeT2",1] <=20)) #  
        # calculate ProbOcc (icloglog(Intercept)) for clearcuts stands:
        ProbOcc[sel,1,] <- icloglog(covdat.col.clearcuts[sel,,1] %*% t(pocc.posteriors.clearcuts) )
    }
    
    ## Mature stands: formula: icloglog(intercept + beta1X1 + beta2X2) = icloglog(pocc.posteriors.mature X covdat.pcol.mature))
    # select based on raw covdat:
    sel <- as.numeric(which(covdat[,"AgeT2",1] >=64)) # standard
    ProbOcc[sel,1,] <- icloglog(covdat.col.mature[sel,,1] %*% t(pocc.posteriors.mature) )
    
    ## Young production stands: 
    # Pocc remains zero, except for c("antrseri", "fomipini", "gloesepi", "tricabie")
    if(myspecies[sp] %in% c("antrseri", "fomipini", "gloesepi", "tricabie")){ # then use OM_mature
        sel <- as.numeric(which(covdat[,"AgeT2",1] > 20 & covdat[,"AgeT2",1] <64  )) # 
        ProbOcc[sel,1,] <- icloglog(covdat.col.mature[sel,,1] %*% t(pocc.posteriors.mature) )
    }
    
    ## also find stands with zero dead wood and set Pocc = 0:
    sel <- as.numeric(which(covdat[,"DWVol_m3perHa_T2_SpruceConif" ,1] ==0)) #
    ProbOcc[sel,1,] <- 0
    

    ## Periods 2-11: simulate based on col-ext models:---------------------
    for(myperiod in 2:11){  # ProbOcc [nstands, nperiods, nsamp]
      
    # make matrices for Pcol and Pext to be filled:
      ep <- matrix(1, nrow = nstands, ncol=nsamp) 
      cp <- matrix(0, nrow = nstands, ncol=nsamp)
      
      ## clearcuts stands: constant probs Pcol = 0, Pext = 1 - except for species with clearcuts stand model
      ## clearcuts stand models for: antrseri, fomipini, gloesepi, phelviti, tricabie
      if(myspecies[sp] %in% c("antrseri", "fomipini", "gloesepi", "phelviti", "tricabie")){ 
        sel <- as.numeric(which(covdat[,"AgeT2",myperiod] <=20))

        # calculate ep, cp from posteriors and covs
        # format: covs[length(sel), ncovs] %*% post[parms, nsamp] = prob[length(sel), nsamp]
        cp[sel,] <- icloglog(covdat.col.clearcuts[sel,,myperiod] %*% t(pcol.posteriors.clearcuts) ) 
        ep[sel,] <- icloglog(covdat.ext.clearcuts[sel,,myperiod] %*% t(pext.posteriors.clearcuts) ) 
      }
      ## for all other species ep=1, cp=0 , as is set up (ie do nothing)
      
      
      ## mature stands: 
      # and YSS if species in: antrseri, fomipini, gloepsei, tricabie
      if(myspecies[sp] %in% c("antrseri", "fomipini", "gloesepi", "tricabie")){ 
        sel <- as.numeric(which(covdat[,"AgeT2",myperiod] >20)) # use mature model also for YSS
      }else{ # for other species set 0
        sel <- as.numeric(which(covdat[,"AgeT2",myperiod] >=64)) # YSS keep ep=1, cp=0
      }
  
      # calculate ep, cp from posteriors and covs
      # format: covs[length(sel), ncovs] %*% post[parms, nsamp] = prob[length(sel), nsamp]
      cp[sel,] <- icloglog(covdat.col.mature[sel,,myperiod] %*% t(pcol.posteriors.mature) ) 
      ep[sel,] <- icloglog(covdat.ext.mature[sel,,myperiod] %*% t(pext.posteriors.mature) ) 
      

      ### calculate occ prob from cp and ep
      ProbOcc[,myperiod,] <- (1 - ProbOcc[,(myperiod-1),])*cp + ProbOcc[,(myperiod-1),]*(1-ep)
      ProbCol[,myperiod,] <- cp   
      ProbExt[,myperiod,] <- ep
       
      # YSS stands: set Pocc to zero: 
      # unless species in: antrseri, fomipini, gloepsei, tricabie
      if(!(myspecies[sp] %in% c("antrseri", "fomipini", "gloesepi", "tricabie"))){ 
         sel <- as.numeric(which(covdat[,"AgeT2",myperiod] > 20  & covdat[,"AgeT2",myperiod] < 64))
         ProbOcc[sel,myperiod,] <- 0
         ProbCol[sel,myperiod,] <- 0
         ProbExt[sel,myperiod,] <- 1
       }
       
      # also find stands with zero dead wood and set Pocc = 0: (else get NaNs from the scaled cov for DWvol=-Inf)
      sel <- as.numeric(which(covdat[,"DWVol_m3perHa_T2_SpruceConif" ,myperiod] ==0)) #
      ProbOcc[sel,myperiod,] <- 0
      
      ## No dead wood: set Pcol to 0 , Pext to 1
      ProbCol[sel,myperiod,] <- 0
      ProbExt[sel,myperiod,] <- 1
       
    }# end periods loop myperiod
    
 
    ## Mechanistic assumption 2: reduce Pocc to 10% in "retention patches" because of edge effects (do this afterwards for all periods)
    sel <- which(standdat$RetentionBufferzone == 1)
    ProbOcc[sel,,] <- 0.1*ProbOcc[sel,,]
    

    ## STORE raw output (optional, large files):
    # use rds file format for single objects
#   saveRDS( ProbOcc, file =paste(OUTdir, myscens[scen],"_", species, "_ProbOccRaw.rds", sep=""))
 
 
    ### Summarize output:-------------------------------------     

    ### Summarise ProbOcc:
    ## Mean etc across all stands: 
    MeanProbOcc_all <- apply(ProbOcc, 2:3, mean)    # sum up across stands to get 11 periods x 1000 samples
    Mean_all <- apply(MeanProbOcc_all, 1, mean) 
    Median_all <- apply(MeanProbOcc_all, 1, median) 
    CRI95_all <- apply(MeanProbOcc_all, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
    rm(MeanProbOcc_all); gc() 
    
    ## Mean etc across production stands: 
    sel <- which(standdat$SetAside == 0)
    MeanProbOcc <- apply( ProbOcc[sel,,] , 2:3, mean)    
    Mean_prod <- apply(MeanProbOcc, 1, mean) 
    Median_prod <- apply(MeanProbOcc, 1, median) 
    CRI95_prod <- apply(MeanProbOcc, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
    rm(MeanProbOcc); gc() 
  
    ## Mean etc across protected stands: 
    sel <- which(standdat$SetAside == 1)
    MeanProbOcc <- apply(ProbOcc[sel,,]  , 2:3, mean) 
    Mean_protect <- apply(MeanProbOcc, 1, mean) 
    Median_protect <- apply(MeanProbOcc, 1, median) 
    CRI95_protect <- apply(MeanProbOcc, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
    rm(MeanProbOcc); gc() 
  
  
    ## collect these in a list:
    results_main <- list(Mean_all= Mean_all, Median_all = Median_all, CRI95_all = CRI95_all,
                       Mean_protect= Mean_protect, Median_protect= Median_protect, CRI95_protect= CRI95_protect,
                       Mean_prod= Mean_prod, Median_prod= Median_prod, CRI95_prod= CRI95_prod)
  
    ## clean up and store summaries:
    rm(Mean_all, Median_all, CRI95_all, Mean_prod, Median_prod, CRI95_prod, Mean_protect, Median_protect, CRI95_protect)
    rm(ProbOcc); gc()
  
    ### Summarise ProbCol:
    MeanProbCol_all <- apply(ProbCol, 2:3, mean)    # sum up across stands to get 11 periods x 1000 samples
    Mean_all <- apply(MeanProbCol_all, 1, mean) 
    Median_all <- apply(MeanProbCol_all, 1, median) 
    CRI95_all <- apply(MeanProbCol_all, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
    rm(MeanProbCol_all); gc() 

    ## Mean etc across production stands: 
    sel <- which(standdat$SetAside == 0)
    MeanProbCol <- apply( ProbCol[sel,,]  , 2:3, mean)   
    Mean_prod <- apply(MeanProbCol, 1, mean) 
    Median_prod <- apply(MeanProbCol, 1, median) 
    CRI95_prod <- apply(MeanProbCol, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
    rm(MeanProbCol); gc() 
  
    ## Mean etc across protected stands: (different weights!)
    sel <- which(standdat$SetAside == 1)
    MeanProbCol <- apply(ProbCol[sel,,]  , 2:3, mean)    # sum up across stands to get 11 periods x 1000 samples
    Mean_protect <- apply(MeanProbCol, 1, mean) 
    Median_protect <- apply(MeanProbCol, 1, median) 
    CRI95_protect <- apply(MeanProbCol, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
    rm(MeanProbCol); gc() 
  
    ## collect these in a list:
    resultsPcol_main <- list(Mean_all= Mean_all, Median_all = Median_all, CRI95_all = CRI95_all,
                       Mean_protect= Mean_protect, Median_protect= Median_protect , CRI95_protect= CRI95_protect,
                       Mean_prod= Mean_prod, Median_prod = Median_prod, CRI95_prod = CRI95_prod)
    rm(Mean_all, Median_all, CRI95_all, Mean_prod, Median_prod, CRI95_prod, Mean_protect, Median_protect, CRI95_protect)
    rm(ProbCol); gc()
  
  
    ### Summarise ProbCol:
    MeanProbExt_all <- apply(ProbExt, 2:3, mean)     
    Mean_all <- apply(MeanProbExt_all, 1, mean) 
    Median_all <- apply(MeanProbExt_all, 1, median) 
    CRI95_all <- apply(MeanProbExt_all, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
    rm(MeanProbExt_all); gc() 
  
    ## Mean etc across production stands:  
    sel <- which(standdat$SetAside == 0)
    MeanProbExt <- apply( ProbExt[sel,,]  , 2:3, mean)   
    Mean_prod <- apply(MeanProbExt, 1, mean) 
    Median_prod <- apply(MeanProbExt, 1, median) 
    CRI95_prod <- apply(MeanProbExt, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
    rm(MeanProbExt); gc() 
  
    ## Mean etc across protected stands: (different weights!)
    sel <- which(standdat$SetAside == 1)
    MeanProbExt <- apply(ProbExt[sel,,]  , 2:3, mean)  
    Mean_protect <- apply(MeanProbExt, 1, mean) 
    Median_protect <- apply(MeanProbExt, 1, median) 
    CRI95_protect <- apply(MeanProbExt, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE) 
    rm(MeanProbExt); gc() 
  
    ## collect these in a list:
    resultsPext_main <- list(Mean_all= Mean_all, Median_all = Median_all, CRI95_all = CRI95_all,
                           Mean_protect= Mean_protect, Median_protect= Median_protect, CRI95_protect= CRI95_protect,
                           Mean_prod= Mean_prod, Median_prod= Median_prod, CRI95_prod= CRI95_prod)
  rm(Mean_all, Median_all, CRI95_all, Mean_prod, Median_prod, CRI95_prod, Mean_protect, Median_protect, CRI95_protect)
  rm(ProbExt); gc()
  
  ## store summary output:--------
  save(results_main, resultsPcol_main,  resultsPext_main, standdat,
       file=paste(OUTdir, myscens[scen],"_", species, "_ProbOcc.RData", sep="")) 
  
  ## clean up:
  rm(results_main, resultsPcol_main,  resultsPext_main) 
  gc()
    
  } # end species loop sp

  
} # end minDiam loop
