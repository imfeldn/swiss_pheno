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
trainyrs <- allyrs[allyrs > 1991]
testyrs <- allyrs[allyrs < 1987]

## early period calibration
# trainyrs <- allyrs[allyrs < 1988]
# testyrs <- allyrs[allyrs > 1991]

print(trainyrs)
print(testyrs)

# --- run calibration for alternating years ----
for(pv in phenovals){

  # read in the full driver files
  # to be subset below to select the
  # correct training and testing years

  # create selection criteria
  # for training and testing
  selection <- data.frame(
    site = stn_list$site,
    year = stn_list$year
    ) |>
    mutate(
      train = year < 1988,
      test = year > 1991
    )

  # select training years
  selection <- selection |>
    dplyr::group_by(site) |>
    mutate(
      train = ifelse(
        length(which(train)) > 1,
        TRUE,
        FALSE
      )
    )

  # select test years
  selection <- selection |>
    dplyr::group_by(site) |>
    mutate(
      train = ifelse(
        length(which(test)) > 1,
        TRUE,
        FALSE
      )
    )

  # subset training and testing datasets
  stn_list_train <- pr_fm_subset(stn_list, selection$train)
  stn_list_test <- pr_fm_subset(stn_list, selection$test)

#   print("fit models")
#
#   mod_out <- pr_fit_comparison(random_seeds = seeds,
#                                models = mods,
#                                data = stn_list_train,
#                                method = "GenSA",
#                                control = list(max.call = 40000, temperature = 1000000000), #/
#                                par_ranges = system.file("extdata", "parameter_ranges.csv", package = "phenor",mustWork = TRUE),
#                                ncores = length(seeds))
#
#   ### set 9999 to NA for evaluation
#   for(ii in 1:length(mod_out$modelled)){
#     mod_out$modelled[[ii]]$predicted_values[mod_out$modelled[[ii]]$predicted_values > 1000] <- NA
#   }
#
#   save(mod_out, file = paste0("data/03_calibrations/",subDir,"/",pv,"/model_comparison_",pv,".RData"))
#   save(stn_list_train, file = paste0("data/03_calibrations/",subDir,"/",pv,"/stn_list_train_",pv,".RData"))
#   save(stn_list_test, file = paste0("data/03_calibrations/",subDir,"/",pv,"/stn_list_test_",pv,".RData"))

}

