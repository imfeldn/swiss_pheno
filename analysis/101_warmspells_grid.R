##### CALCULATION OF WARM SPELL DURATION INDEX #####

#### packages
library(ncdf4);library(stats);library(geosphere);library(tictoc);library(lubridate)
library(foreach);library(iterators);library(parallel);library(doParallel)

source("R/helpfuns_temp_quantiles.R")

#### directories
indir <- "../swiss_recon/temp/"
outdir <- "data/01_climate_indices/wsdi/"

### temperature data
tempfile <- list.files("../swiss_recon/temp/", pattern = "CH_temp", full.names = T)
ft <- nc_open(tempfile[1])

### select quantile
qth <- 0.9

### calculate for every grid cell and every day of the year the selcted quantile
calc <- F

if(calc){

  ## load file with pre-selected reference period
  nc <- nc_open(paste0(indir,"clim/CH_temp_EnKF_1871-01-01-1900-12-31.nc"))
  time <- as.Date(ncvar_get(nc,"time"), origin= "1971-01-01")
  lon <- ncvar_get(nc,"E")
  lat <- ncvar_get(nc,"N")
  qtemp <- array(NA, dim = c(length(lon),length(lat),365))

  for(ll in 1:length(lon)){
    print(ll)
    dat <- ncvar_get(nc,"temp", start = c(ll,1,1), count = c(1,-1,-1))
    qtemp[ll,,] <- calc_quants_grid(temp = dat, date = time, qth = qth)
  }
  save(qtemp, file = paste0("data/01_climate_indices/wsdi/qtemp_th0",qth*10,"_1871-1900.RData"))
} else{

  load(paste0("data/01_climate_indices/wsdi/qtemp_th0",qth*10,"_1871-1900.RData"))
}

### function for warm spells
warmspellindex <- function(temp, dates, qdat){

  if(all(is.na(temp))){
    spelldays <- rep(NA,4)
  }
  else{
    jdays <- as.numeric(as.character(get_date_factors_obj_special(dates)$jdays))
    date_fac <- as.factor(ifelse(month(dates) %in% c(3:5),"MAM",
                                  ifelse(month(dates) %in% c(6:8),"JJA",
                                         ifelse(month(dates) %in% c(9:11),"SON","DJF"))))

    f <- match.fun(">")
    min_length <- 6
    d <- f(temp, qdat[jdays])
    ## set to true if only one day in between is FALSE
    d[which(c(FALSE,d[1:length(d)-1]) & !d[1:length(d)] & c(d[2:length(d)],FALSE))] <- TRUE

    d2 <- Reduce(function(x, y) { # Reduce sucessively combines elements of a given vector
      z <- c(rep(FALSE, y), d[1:(length(d) - y)]) & x
      return(z)
    }, 1:(min_length-1), d)

    periods <- Reduce(function(x, y) {
      return(c(d2[(y + 1):length(d2)], rep(FALSE, y)) | x)
    }, 1:(min_length-1), d2)

    spelldays  <- climdex.pcic:::tapply.fast(periods, date_fac, sum)

    ## set to NA, if less than 3 month of data
    fac_len <- climdex.pcic:::tapply.fast(date_fac,date_fac,length)
    spelldays[fac_len < 85] <- NA

  }
  return(spelldays)
}


###  create new files
print("create file")
filename1 <- paste(outdir,"/WSDI_season_1763-2020_",Sys.Date(),".nc",sep="")

N <- ncvar_get(ft, "N")
E <- ncvar_get(ft, "E")
nx <- length(E)
ny <- length(N)
lon1 <- ncdim_def("E", "meter", E)
lat2 <- ncdim_def("N", "meter", N)

dates <- seq(as.Date("1763-01-01"),as.Date("2020-12-31"), by = "day")
time_hist <- which(substr(dates,6,10) %in% c("01-01","04-01","07-01","10-01")) - which(dates=="1970-01-01")
time_def <- ncdim_def("time","days since 1970-01-01", time_hist, unlim=TRUE)
mv <- -999 #missing value to use
var_temp1 <- ncvar_def("days", "days", list(lon1, lat2, time_def), longname="ndays heat-wave", mv)
ncnew1 <- nc_create(filename1, list(var_temp1))

### register cores
registerDoParallel(cores = floor(detectCores()*0.40))
ff <- which(grepl("1961",tempfile))

### loop over years
print("start calc")
for(ff in 1:length(tempfile)){
  ft <- nc_open(tempfile[ff])
  ty <- ncvar_get(ft,varid=switch(grepl("EnKF",tempfile[ff]) + 1,"TabsD","temp"))
  time <- as.Date(ncvar_get(ft,varid="time"), origin = switch(grepl("EnKF",tempfile[ff]) + 1,"1900-01-01","1970-01-01") )
  nc_close(ft)

  if(ff != length(tempfile)){
    ft2 <- nc_open(tempfile[ff+1])
    ty2 <- ncvar_get(ft2,varid=switch(grepl("EnKF",tempfile[ff+1]) + 1,"TabsD","temp"))
    time2 <- as.Date(ncvar_get(ft2,varid="time"),origin = switch(grepl("EnKF",tempfile[ff]) + 1,"1900-01-01","1970-01-01"))
    nc_close(ft2)

    ## merge to be able to calculate seasonal statistics
    newtime <- c(time,time2)
    mars <- which(substr(newtime,6,10) == "03-01")
    newdat <- abind::abind(ty,ty2)[,,mars[1]:(mars[2]-1)]
    newtime <- newtime[mars[1]:(mars[2]-1)]

  } else {
    newtime <- time
    mars <- which(substr(newtime,6,10) == "03-01")
    newdat <- ty[,,mars[1]:length(time)]
    newtime <- newtime[mars[1]:length(time)]
  }

  temp_list <- split(newdat,seq(nrow(newdat)*ncol(newdat)))
  qtemp_list <- split(qtemp,seq(nrow(qtemp)*ncol(qtemp)))

  print(paste0("recon ",ff + 1762))
  tictoc::tic()
  out_list <- foreach(a = 1:length(temp_list)) %dopar% {
    warmspellindex(temp = temp_list[[a]], dates = newtime, qdat = qtemp_list[[a]])
  }
  tictoc::toc()

  spelldays <- aperm(array(unlist(out_list), dim=c(4,nrow(newdat),ncol(newdat))), perm = c(2,3,1))

  ### change order of seasons
  spelldays <- spelldays[,,c(3,2,4,1)]
  ff2 <- ((year(time)[1] - 1763 + 1) * 4) - 2 ## select spring of the year

  ### write to file
  ncvar_put(ncnew1, var_temp1, spelldays,  start = c(1,1,ff2), count = c(-1,-1,4))
 }

nc_close(ncnew1)
print("done")
