#### CALCULATE INDICES FOR 20CRV3 ########

### packages
library(ncdf4);library(lubridate);library(mmand);library(doParallel)

### helpfer functions
source("R/helpfuns.R")
source("R/helpfuns_swissplateau.R")
source("R/helpfuns_temp_quantiles.R")

### function to restructure 20cr
func <- function(dat, ...){
  args <- list(...)

  if(length(args$qdat) > 0) {
    out <- sapply(1:nrow(dat), function(x){do.call(get(indices[ii]), c(list(dat[x,],args$qdat[x,],args)))})
  } else{
    out <- sapply(1:nrow(dat), function(x){do.call(get(indices[ii]), c(list(dat[x,],args)))})
  }

  yrs <- length(unique(year(args$date)))

  if(nrow(out) > (yrs * 2)){
    ll <- length((yrs+1):nrow(out))
    dims <- c(yrs,ll/yrs,dt[1],dt[2])
    perms <- c(3,4,1,2)
  } else {
    dims <- c(length((yrs + 1):nrow(out)),dt[1],dt[2])
    perms <- c(2,3,1)
  }

  out2 <- aperm(array(out[(yrs+1):nrow(out),], dim = dims), perm = perms)
  time <- out[1:yrs,1]
  return(list(ind = out2,time = time))
}

### indices
indices <- c("gsstart","gsend","gsl","ffd","lfd","fd","gdd","summerdays","wsdi","csdi","tempmean")

### read data
path <- "../swiss_indices/2_data/20crv3/t2m/swiss_t2m_daily_1806-2015_new.nc"
nc <- nc_open(path)
temp <- ncvar_get(nc, "t2m") - 273.15
lon <-  ncvar_get(nc, "longitude")
lat <-  ncvar_get(nc, "latitude")
time <- as.Date(ncvar_get(nc, "time")/24, origin = "1800-01-01")

## resturcture temp
dt <- dim(temp)
temp_mat <- array(temp, dim = c(dt[1] * dt[2], dt[3]))

index_list <- list()
ii <- 6
thval <- 200

registerDoParallel(cores=floor(detectCores()*0.35))

for(ii in 1:length(indices)){

  print(indices[ii])

  if(indices[ii] == "gdd"){
    index_list[[ii]] <- func(dat = temp_mat, date = time, th = thval)
  } else if(indices[ii] %in% c("csdi","wsdi")){
    qth <- ifelse(indices[ii] == "csdi", 0.1, 0.9)
    qdat <- calc_quants_grid(temp = temp_mat, date = time, qth = qth)
    index_list[[ii]] <- func(dat = temp_mat, dates = time, qdat = qdat)
  } else {
    index_list[[ii]] <- func(dat = temp_mat, date = time)
  }
}

### add precipitation indices as comparison

names(index_list) <- indices
index_list_20crv3 <- list(index=index_list,lon=lon,lat=lat)
save(index_list_20crv3, file = paste0("data/01_climate_indices/20crv3_indices_",Sys.Date(),".RData"))

