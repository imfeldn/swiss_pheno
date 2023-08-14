############# RECONSTRUCT PAST PHENOLOGY BASED ON ONE MODEL CALIBRATION #########
rm(list=ls())

#### packages
library(ncdf4);library(geosphere);library(tictoc)
library(foreach);library(phenor)

#### directories
outdir <- "data/02_pheno_recon/"

### load best models
annealing_version <- "sa_run_2023-08-09_class123"

#pv <- "maesh13d"

for (pv in phenovals){

  print(pv)

  (load(paste0("data/03_calibrations/",annealing_version,"/",pv,"/aicc_mods_",pv,".RData")))

  modsel <- which(aicc_mat == min(aicc_mat), arr.ind = T)
  (load(list.files(paste0("data/03_calibrations/",annealing_version,"/",pv),pattern = "model_comparison",full.names = T)))
  opt_params <- mod_out$modelled[[modsel[1]]]$parameters[modsel[2],]

  ### temperature data
  tempfile <- list.files("../swiss_recon/temp/", pattern = "CH_temp_", full.names = T)
  ft <- nc_open(tempfile[1])

  ##  create new files
  filename <- paste(outdir,"/CH_",rownames(modsel),"_",pv,"_1763-01-01-2020-12-31_",Sys.Date(),".nc",sep="")

  N <- ncvar_get(ft, "N");E <- ncvar_get(ft, "E")
  nx <- length(E);ny <- length(N)
  lon1 <- ncdim_def("E", "meter", E);lat2 <- ncdim_def("N", "meter", N)

  ### get lat values for estimating the daylength
  latdegb <- matrix(NA,nrow=length(E),ncol=length(N))
  Ep <- (E-2600000)/1000000
  Np <- (N-1200000)/1000000

  for (i in 1:length(E)){
    for(j in 1:length(N)){
    latdegb[i,j] <- (16.9023892+3.238272*Np[j]-0.270978*(Ep[i]^2)-0.002528*(Np[j]^2)-0.0447*(Ep[i]^2)*Np[j]-0.0140*(Np[j]^3))*100/36
  }}

  latvec <- as.vector(latdegb)
  Lidoy <- sapply(latvec, function(x) phenor::daylength(1:365,x))

  dates <- seq(as.Date("1763-01-01"),as.Date("2020-12-31"), by = "day")
  time_hist <- which(substr(dates,6,10) == "01-01")
  time_def <- ncdim_def("time","days since 1763-01-01", time_hist, unlim=TRUE)
  mv <- -999 #missing value to use
  var_temp <- ncvar_def("DOY", "doy", list(lon1, lat2, time_def), longname=pv, mv)
  ncnew <- nc_create(filename, list(var_temp),force_v4=T)
  ncatt_put(nc = ncnew, varid = "DOY", attname = "model", attval = modsel, prec=NA)
  ncatt_put(nc = ncnew, varid = "DOY", attname = "params", attval = paste0(opt_params,collapse="|"), prec=NA)

  func <- get(rownames(modsel))

  registerDoParallel(cores = floor(detectCores()*0.45))

  for(ff in 1:(length(tempfile)-1)){

    ft1 <- nc_open(tempfile[ff])
    time1 <- as.Date(ncvar_get(ft1, "time"), origin = switch(grepl("EnKF",tempfile[ff]) + 1,"1900-01-01","1970-01-01"))
    ty1 <- ncvar_get(ft1,varid = switch(grepl("EnKF",tempfile[ff]) + 1,"TabsD","temp"))[,,(length(time1) - 101):length(time1)]
    nc_close(ft1)

    ft2 <- nc_open(tempfile[ff + 1])
    ty2 <- ncvar_get(ft2,varid = switch(grepl("EnKF",tempfile[ff + 1]) + 1,"TabsD","temp"))[,,1:263]
    nc_close(ft2)

    tynew <- abind::abind(ty1,ty2)
    temp_list <- split(tynew,seq(nrow(tynew) *ncol(tynew)))

    doy <- c(-102:-1,1:263)
    temp_list <- lapply(1:length(temp_list), function(x) y <- list(Ti = matrix(temp_list[[x]], ncol = 1), Li = matrix(Lidoy[c(264:365,1:263),x], ncol = 1), doy = doy))

    print(paste0("recon ",ff + 1763))
    tictoc::tic()
    doy_list <- foreach(a = 1:length(temp_list)) %dopar% {
      func(par = opt_params, data = temp_list[[a]])
    }
    tictoc::toc()

    doymat <- matrix(unlist(doy_list), nrow = nrow(tynew), ncol = ncol(tynew))
    doymat[doymat == 9999] <- NA
    ncvar_put(ncnew, var_temp, doymat,  start = c(1,1,ff + 1), count = c(-1,-1,1))
  }
  nc_close(ncnew)
  print("done")
}
