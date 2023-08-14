############# CALCULATE FROST INDEX FOR DIFFERENT PHENOPHASES ################
rm(list=ls())

### packages
library(ncdf4);library(parallel);library(foreach)
library(abind);library(doParallel)

### directories
plotdir <- "data/04_pheno_recon/"

### data path
tempfiles <- list.files("../swiss_recon/temp/", pattern = "CH_temp_", full.names = T)

### phenology file
phenofiles <- list(mprua65d = "data/04_pheno_recon/CH_PTTs_mprua65d_1763-01-01-2020-12-31_2023-08-02.nc")


th <- 0
pv <- 1

## read pheno
for(pv in 1:length(phenofiles)){

  print(names(phenofiles)[pv])
  nc <- nc_open(phenofiles[[pv]])
  pheno <- ncvar_get(nc,"DOY")
  pheno_N <- ncvar_get(nc,"N")
  pheno_E <- ncvar_get(nc,"E")
  pheno[pheno == 1] <- NA
  pheno[,,1] <- NA

  newpheno <- pheno - 5 ### chose 5 days before phenology

  ## create nc files with number of days after pheno blossom
  tmean_after_pheno <- array(NA, dim(pheno))

  cores <- 5
  registerDoParallel(cores=cores)
  i <- 20

  for(i in 1:length(tempfiles)){
    print(yr <- i + 1762)
    nc <- nc_open(tempfiles[i])
    tval <- ncvar_get(nc,ifelse(yr < 1961,"temp","TabsD"))[,,1:212]
    temp_0 <- aperm(apply((tval <= 0) + 0 ,1:2,"*", c(1:212) ), c(2,3,1)) ## this turns array

    # sum of temperature with tmean below 0 °C 30 days after flowering
    temp_af <- mclapply(1:nrow(tval), function(x) {
      temp2d = tval[x,,]
      temp2d[(temp_0[x,,] < newpheno[x,,i])] = NA ## go five days back to include frost events slightly prior to the phenophase as well
      ll = sapply(1:nrow(temp2d),function(j) if(!all(is.na(temp2d[j,]))){
        if(!is.na(pheno[x,j,i])){
          cherrind <- newpheno[x,j,i]:212 # select all doys from onset to end of July
          td <- temp2d[j,cherrind]
          sum(td[td <= th], na.rm = T)
        } else {NA}
      } else {NA})
      return(ll)
    }, mc.cores = cores)
    temp_af <- do.call(rbind, temp_af)
    tmean_after_pheno[,,i] <- temp_af

    nc_close(nc)
  }

  ### wrtie files to folder
  save( tmean_after_pheno,
       file = paste0("data/04_pheno_recon/CH_frostindex_th",th,"_after_",names(phenofiles)[pv],"_5daysprior_",Sys.Date(),".RData"))

  print("files saved")

}
