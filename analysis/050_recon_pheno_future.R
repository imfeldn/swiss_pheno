rm(list=ls())

setwd("/scratch3/noemi/wear/swiss_indices/")

library(ncdf4);library(doParallel);library(lubridate);library(abind);library(geosphere)

path <- "/mnt/climstor/CH2018/QMgrid/tas/"
outdir <- "/scratch3/noemi/wear/swiss_indices/5_output/5_pheno_future/"

files <- list.files(path,pattern = "CH2018", full.names = T)

swiss_mean_latitude <- 46.795953732138
dl <-  daylength(swiss_mean_latitude,1:360)

mod <- "gddpt1"
pv <- "mfags13d"

## params
(load(paste0("5_output/3_calibrations/1_annealing/",mod,"_optim_param_",pv,"_2022-12-22.RData")))

#### function ####

## interpolate the models with only 360 days
intrp <- function(temp){
  missdays <- which(apply(temp,3, function(x){all(is.na(x))}))
  missdays <- missdays[missdays < 364] ## just interpolate first few days 
  for(ii in missdays){
    if(ii != 1){
      temp[,,ii] <- (temp[,,ii-1] + temp[,,ii+1])/2
    } else {
      temp[,,ii] <- temp[,,ii+1] ### this is a work-around if the 1jan is missing - if it's decided what to do, this should be changed to something better.
    }
    }
  return(temp)
}

## get date
gddpt1_fun <- function(temp, params){
    Fstar = params[1]; dt = params[2]; Rp = params[3]
    ## interpolate the missing values
    temp_full <- intrp(temp)
    gd <- apply(temp_full[,,1:length(dl)],1:3,"-",dt)
    gd[gd<0] <- 0
    fd <- 1-(0.002*Rp*((20-dl)^2))
    gp <- apply(gd,1:2,"*",fd)
    gp[gp<0] <- 0
    gs1 <- apply(gp,2:3,cumsum)
    gs2 <- gs1 > Fstar
    gdd2 <- apply(gs2,2:3,function(x) {
      if(any(is.na(x))){
        dd = NA ## set models with missing dates to NA instead of 1
        } else {
        dd = which.max(x)
        if(length(dd) == 0){dd = NA} ## set grid cells which do not reach threshold to NA
      }
      return(dd)
      })
    return(gdd2)
}

## get date
sigmoidpt1_fun <- function(temp, params){
  Fstar =  params[1]; T50 =  params[2]; dt =  params[3]; Rp = params[4]
  ## interpolate the missing values
  temp_full <- intrp(temp)
  rf_help <- apply(temp_full[,,1:length(dl)],1:3,function(y) {1 + (exp(-dt*(y-T50)))})
  rf <- 1/rf_help
  rf[temp_full[,,1:length(dl)] < 0] <- 0 ## set values to zero, that are below zero degrees
  
  fd <- c(1-(0.002*Rp*((20-dl)^2)))[1:length(dl)]
  rf <- apply(rf, 1:2,"*",fd)
  gs1 <- apply(rf,2:3,cumsum)
  gs2 <- gs1>Fstar
  
  gdd2 <- apply(gs2,2:3,function(x) {
    if(any(is.na(x))){
      dd = NA ## set models with missing dates to NA instead of 1
    } else {
      dd = which.max(x)
      if(length(dd) == 0){dd = NA} ## set grid cells which do not reach threshold to NA
    }
    return(dd)
  })
  return(gdd2)
}

func <- get(paste0(mod,"_fun"))
registerDoParallel(cores=floor(detectCores()*0.16)) 
acomb <- function(...) abind(..., along=3)


ff <- 34
#files <- files[nanmodels]

for(ff in 1:length(files)){

  print(files[ff])
  
  nc <- nc_open(files[ff])
  lat <- ncvar_get(nc, "lat")
  lon <- ncvar_get(nc, "lon")
  time <- as.Date(ncvar_get(nc,"time"), origin = "1900-01-01")
  
  ch_mod <- gsub(".nc","",gsub(".*CH2018_tas_","",files[ff]))
  ## create new file for beech leaf unfolding CH2018
  filename <- paste(outdir,"/CH2018_",pv,"_",ch_mod,"_",Sys.Date(),".nc",sep="")
  
  nx <- length(lon)
  ny <- length(lat)
  lon1 <- ncdim_def("lon", "degrees_east", lon)
  lat2 <- ncdim_def("lat", "degrees_north", lat)
  
  time_sel <- which(substr(time,6,10) == "01-01")
  time_def <- ncdim_def("time","days since 1981-01-01", time_sel, unlim=TRUE)
  mv <- -999 #missing value to use
  var_temp <- ncvar_def("BLU", "day of year", list(lon1, lat2, time_def), longname="Beech leaf unfolding", mv) 
  ncnew <- nc_create(filename, list(var_temp))

  doy <- array(NA, dim = c(nx,ny,length(time_sel)))
  pers <- round(seq(1,length(time_sel),length.out = 10))
  
  print("recon leaf unfolding")
  
  for(p in 1:(length(pers)-1)){
    
    ifelse(p == (length(pers)-1),
           tind <- time_sel[pers[p]]:length(time),
           tind <- time_sel[pers[p]]:(time_sel[pers[p+1]]-1))
    ifelse(p == (length(pers)-1),
           tind2 <- pers[p]:(pers[p+1]),
           tind2 <- pers[p]:(pers[p+1]-1))
    tas <- ncvar_get(nc = nc, varid = "tas", start = c(1,1,tind[1]), count = c(-1,-1,length(tind)))
    time2 <- as.Date(ncvar_get(nc = nc, varid = "time", start = c(tind[1]), count = c(length(tind))), origin = "1900-01-01")
    
    tictoc::tic()
    doy[,,tind2] <- foreach(a = 1:length(tind2), .combine = "acomb") %dopar% {
      aa = year(time2) %in% year(time[time_sel[tind2[a]]])
      bb = func(temp = tas[,,aa], params = opt_params[2:length(opt_params)])
      return(bb)  
    }
    tictoc::toc()
  }
  
  print("write leaf unfolding")
  ncvar_put(ncnew, var_temp, doy)
  nc_close(ncnew)
  
}

print("done")
