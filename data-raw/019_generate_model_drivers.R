# Generate model drivers from climate
# reconstruction data, will be subset
# during model calibration into a
# test and training datasets

# load packages
library(phenor)

# set years to cover
allyrs <- 1951:2020

# select phenophases for calibration
# phenovals <- c("mfags13d","mprua65d","maesh13d",
#                "mtaro65d","manen65d","mmald65d",
#                "mpyrc65d","mcora65d","mcora13d",
#                "mtusf65d","mcarp65d","mlard13d")

phenovals <- c("maesh13d","mfags13d","mprua65d")

# where to save the data
output_dir <- "data/02_pheno_net/"

# create the directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Generating complete driver data (all years, no subsets)")
# --- run calibration for alternating years ----
for(pv in phenovals){

  message(paste0("-- processing: ", pv))

  # set filenames for driver data
  meteo_file <- file.path("data-raw/02_pheno_net/",
                          paste0("meteo_",pv,"_class123.RData")
                          )
  phenology_file <- file.path("data-raw/02_pheno_net/",
                           paste0("meteo_",pv,"_class123.RData")
                           )

  # load pheno data and meteorological data
  load(meteo_file)
  load(phenology_file)

  # split out station names
  stns_pv <- unique(pheno_meta$nat_abbr)

  # create empty list to hold all data
  stn_list <- list()

  # convert data to phenor format
  for(stn in 1:length(stns_pv)){

    # get id to process
    id <- stns_pv[stn]

    ind <- which(pheno_meta$nat_abbr == id & pheno_meta$reference_year < 2021)
    ltm <- apply(meteo_dat[rownames(meteo_dat) == id,,], 2, mean)

    # create meteo data subset
    meteo_subset <- meteo_dat[rownames(meteo_dat) == id, allyrs %in% pheno_meta$reference_year[ind],]

    # check if this is a matrix
    # if not subsetting the array
    # failed and the site does not exist
    # (and should be skipped)
    if(nrow(meteo_subset) == 0){
      stn_list[[stn]] <- NULL
    } else {

    locs <- as.numeric(unlist(pheno_meta[ind[1],c("lat","lon")]))

    stn_list[[stn]]  <- list(
      site = id,
      location = locs,
      doy = c(-102:-1,1:263),
      ltm = ltm, # long term mean temperature for a given location / ignored
      transition_dates = pheno_meta$doy[ind],
      year = pheno_meta$reference_year[ind],
      Ti = t(t_mean),
      Tmini = NULL, # tmin and max ignored / allowed to be empty
      Tmaxi = NULL, #
      Li = matrix(
        rep(daylength(c(264:365,1:263),locs[1]),
            length(ind)),
        nrow = 365,
        ncol = length(ind)
        )
      )
    }
  }

  # add names to the list elements
  names(stn_list) <- stns_pv

  # flatten the phenocam nested list format
  # for easier subsetting
  stn_list <- pr_flatten(stn_list)

  # get class (based on meteo data)
  basename(meteo_file)

  # construct filename
  filename <- file.path(output_dir, paste0(pv, "_full_model_drivers.rds"))

  # save data as a compressed RDS file
  saveRDS(
    stn_list,
    filename,
    compress = "xz"
    )
}
