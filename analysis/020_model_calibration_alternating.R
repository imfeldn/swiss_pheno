# Calibration routine using alternating years
# over the whole dataset
try(detach("package:phenor", unload = TRUE))
### load packages
library(phenor)
library(dplyr)

### select phenophases for calibration
# phenovals <- c("mfags13d","mprua65d","maesh13d",
#                "mtaro65d","manen65d","mmald65d",
#                "mpyrc65d","mcora65d","mcora13d",
#                "mtusf65d","mcarp65d","mlard13d")

phenovals <- c("maesh13d","mfags13d","mprua65d")
classes <- "123_longseries"

### create directory for new calibration with simulated annealing
base_dir <- file.path(
  "data/03_calibrations/",
  paste0("sa_run_",Sys.Date(),"_class",classes)
)

dir.create(
  base_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

### select seeds, model, years
#seeds <- c(1,7,22,31,45,57,63,99)
seeds <- 1
mods <-  c("LIN", "TT", "TTs", "PTT", "PTTs","M1","M1s","AT")

# list driver files
driver_files <- list.files("data/02_pheno_net/","*.rds", full.names = TRUE)

message("Calibrating models:")
# --- run calibration for alternating years ----
lapply(phenovals, function(pv) {

  message(paste0("-- calibrating: ", pv))

  # read in the full driver files
  # to be subset below to select the
  # correct training and testing years
  drivers <- readRDS(driver_files[grepl(paste0(pv,".*",classes,".rds"), driver_files)])

  # create selection criteria
  # for training and testing
  selection <- data.frame(
    site = drivers$site,
    year = drivers$year
  ) |>
    mutate(
      train = (year %% 2 == 0),
      test = (year %% 2 == 1)
    )

  # select training years
  selection <- selection |>
    dplyr::group_by(site) |>
    mutate(
      train = ifelse(
        length(which(train)) > 1 & train == TRUE,
        TRUE,
        FALSE
      )
    )

  # select test years
  selection <- selection |>
    dplyr::group_by(site) |>
    mutate(
      test = ifelse(
        length(which(test)) > 1 & test == TRUE,
        TRUE,
        FALSE
      )
    )

  # subset training and testing datasets
  stn_list_train <<- pr_fm_subset(drivers, selection$train)
  stn_list_test <- pr_fm_subset(drivers, selection$test)

  print("fit models")
  nr_seeds <- length(seeds)

  mod_out <- pr_fit_comparison(
    random_seeds = 1,
    models = mods,
    data = stn_list_train,
    method = "bayesiantools",
    control = list(
      sampler = "DEzs",
      settings = list(
        burnin = 1000,
        iterations = 4000
        )
    )
  )

  ### set 9999 to NA for evaluation
  for(ii in 1:length(mod_out$modelled)){
    mod_out$modelled[[ii]]$predicted_values[mod_out$modelled[[ii]]$predicted_values > 1000] <- NA
  }

  # create output directory
  dir.create(
    file.path(base_dir, pv),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # save data as compressed RDS files
  saveRDS(
    mod_out,
    file = file.path(base_dir,pv,paste0("model_comparison_",pv,".rds")),
    compress = "xz"
  )

  saveRDS(
    stn_list_train,
    file = file.path(base_dir,pv,paste0("stn_list_train_",pv,".rds")),
    compress = "xz"
  )

  saveRDS(
    stn_list_test,
    file = file.path(base_dir,pv,paste0("stn_list_test_",pv,".rds")),
    compress = "xz"
  )
})
