##### CALIBRATIONS OF DIFFERENT PROCESS-BASED PHENOLOGICAL MODELS #######

### load packages
library(phenor)

### select phenophases for calibration
# phenovals <- c("mfags13d","mprua65d","maesh13d",
#                "mtaro65d","manen65d","mmald65d",
#                "mpyrc65d","mcora65d","mcora13d",
#                "mtusf65d","mcarp65d","mlard13d")

phenovals <- c("maesh13d","mfags13d","mprua65d")

### create directory for new calibration with simulated annealing
subDir <- paste0("sa_run_",Sys.Date(),"_class123_calibration_on_earlyperiod")
dir.create(file.path("data/03_calibrations/",subDir), showWarnings = F)

### select seeds, model, years
seeds <- c(1,7,22,31,45,57,63,99)
mods <-  c("LIN", "TT", "TTs", "PTT", "PTTs","M1","M1s","AT")
allyrs <- 1951:2020

trainyrs <- allyrs[allyrs %% 2 == 0]
testyrs <- allyrs[allyrs %% 2 == 1]

## late period calibration
#trainyrs <- allyrs[allyrs > 1991]
#testyrs <- allyrs[allyrs < 1987]

## early period calibration
trainyrs <- allyrs[allyrs < 1988]
testyrs <- allyrs[allyrs > 1991]

print(trainyrs)
print(testyrs)

### run calibration for all phenophases
for(pv in phenovals){

  print(pv)

  # create subdir for pv
  dir.create(file.path("data/03_calibrations/",subDir,"/",pv), showWarnings = F)

  # load pheno data and meteorological data
  (load(paste0("data/02_pheno_net/meteo_",pv,"_class123.RData")))
  (load(paste0("data/02_pheno_net/pheno_meta_",pv,"_class123.RData")))

  stns_pv <- unique(pheno_meta$nat_abbr)
  stn_list <- list()

  # convert data to phenor format
  for(stn in 1:length(stns_pv)){

    id <- stns_pv[stn]
    ind <- which(pheno_meta$nat_abbr == id & pheno_meta$reference_year < 2021)
    locs <- as.numeric(unlist(pheno_meta[ind[1],c("lat","lon")]))
    stn_list[[stn]]  <- list(site = id, location = locs,  doy = c(-102:-1,1:263),
                             ltm = rep(15,365), # I'm not sure what this is for
                             transition_dates = pheno_meta$doy[ind],
                             year = pheno_meta$reference_year[ind],
                             Ti = t(meteo_dat[rownames(meteo_dat) == id,allyrs %in% pheno_meta$reference_year[ind],]),
                             Tmini = t(meteo_dat[rownames(meteo_dat) == id,allyrs %in% pheno_meta$reference_year[ind],]) - 3, # this is artificial data for tmin, but never used in the models
                             Tmaxi = t(meteo_dat[rownames(meteo_dat) == id,allyrs %in% pheno_meta$reference_year[ind],]) + 3, # this is artificial data for tmax, but never used in the models
                             Li = matrix(rep(daylength(c(264:365,1:263),locs[1]), length(ind)), nrow = 365, ncol = length(ind)))


  }

  names(stn_list) <- stns_pv

  ### select training and testing data
  stn_list_train <- lapply(stn_list, function(x){
    z <- x
    selyrs <- x$year%in% trainyrs
    z$transition_dates <- x$transition_dates[selyrs]
    z$year <- x$year[selyrs]
    z$Ti <- x$Ti[,selyrs]
    z$Tmini <- x$Tmini[,selyrs]
    z$Tmaxi <- x$Tmaxi[,selyrs]
    z$Li <- x$Li[,selyrs]
    return(z)
  })

  ### only select time series with more than 1 year of data
  stn_list_train <- stn_list_train[sapply(stn_list_train, function(x) length(x$year)) > 1]

  stn_list_test <- lapply(stn_list, function(x){
    z <- x
    selyrs <- x$year %in% testyrs
    z$transition_dates <- x$transition_dates[selyrs]
    z$year <- x$year[selyrs]
    z$Ti <- x$Ti[,selyrs]
    z$Tmini <- x$Tmini[,selyrs]
    z$Tmaxi <- x$Tmaxi[,selyrs]
    z$Li <- x$Li[,selyrs]
    return(z)
  })

  stn_list_test <- stn_list_test[sapply(stn_list_test, function(x) length(x$year)) > 1]

  print("fit models")

  mod_out <- pr_fit_comparison(random_seeds = seeds,
                               models = mods,
                               data = stn_list_train,
                               method = "GenSA",
                               control = list(max.call = 40000, temperature = 1000000000), #/
                               par_ranges = system.file("extdata", "parameter_ranges.csv", package = "phenor",mustWork = TRUE),
                               ncores = length(seeds))

  ### set 9999 to NA for evaluation
  for(ii in 1:length(mod_out$modelled)){
    mod_out$modelled[[ii]]$predicted_values[mod_out$modelled[[ii]]$predicted_values > 1000] <- NA
  }

  save(mod_out, file = paste0("data/03_calibrations/",subDir,"/",pv,"/model_comparison_",pv,".RData"))
  save(stn_list_train, file = paste0("data/03_calibrations/",subDir,"/",pv,"/stn_list_train_",pv,".RData"))
  save(stn_list_test, file = paste0("data/03_calibrations/",subDir,"/",pv,"/stn_list_test_",pv,".RData"))

}

