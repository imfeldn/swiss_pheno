# Generate model drivers from climate
# reconstruction data, will be subset
# during model calibration into a
# test and training datasets

# load packages
library(phenor)

# select phenophases for calibration
# phenovals <- c("mfags13d","mprua65d","maesh13d",
#                "mtaro65d","manen65d","mmald65d",
#                "mpyrc65d","mcora65d","mcora13d",
#                "mtusf65d","mcarp65d","mlard13d")

phenovals <- c("maesh13d","mfags13d","mprua65d")

output_dir <- file.path("data/02_pheno_net/",subDir)
dir.create(output_dir, showWarnings = F)

# --- run calibration for alternating years ----
for(pv in phenovals){

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

  # add names to the list elements
  names(stn_list) <- stns_pv

  # flatten the phenocam nested list format
  # for easier subsetting
  stn_list <- pr_flatten(stn_list)

  # construct filename
  filename <- file.path(output_dir, paste0(pv, "_full_model_drivers.rds"))

  # save data as a compressed RDS file
  saveRDS(
    stn_list,
    filename,
    compress = "xz"
    )
}

