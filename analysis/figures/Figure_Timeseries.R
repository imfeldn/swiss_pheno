
rm(list=ls())

setwd("/scratch3/noemi/wear/swiss_indices/")
library(ncdf4);library(lubridate);library(mmand)

source("/scratch3/noemi/ncode/R/gaussian_wrap.R")
source("1_code/R/000_helpfuns.R")
plotdir <- "5_output/2_plots/1_climate_indices/all/"

(load("5_output/1_climate_indices/swissplateau_stations_indices_2023-07-14.RData"))
(load("5_output/1_climate_indices/20crv3_indices_2023-07-28.RData"))

### load wet days
nc <- nc_open("2_data/20crv3/prcptot/WD_season_1806-2015.nc")
wds_20cr <- ncvar_get(nc,"tp")
lon_wds <-  ncvar_get(nc,"longitude");lat_wds <-  ncvar_get(nc,"latitude")
time_wds <- as.Date(ncvar_get(nc,"time")/24, origin = "1800-01-01")

swisslats <- c(46.1,47.8)
swisslons <- c(5.7,6.9)

iind <- which(index_list_20crv3$lon >= swisslons[1] & index_list_20crv3$lon <= swisslons[2])
jind <- which(index_list_20crv3$lat  >= swisslats[1] & index_list_20crv3$lat <= swisslats[2])[1]

iind_wd <- which(lon_wds >= swisslons[1] & lon_wds <= swisslons[2])
jind_wd <- which(lat_wds  >= swisslats[1] & lat_wds <= swisslats[2])[1]

#iind <- 7;jind <- 8
index_list_20crv3$lon[iind]
index_list_20crv3$lat[jind]

sfc <- c("Alp", "AlpsS","Jura", "Swiss Plateau", "Pre-Alps")
regind <- which(sfc == "Swiss Plateau")

indices <- list(tempmean = list(file = "5_output/1_climate_indices/tempmean/fldmean_regions_TEMPMEAN_season_1763-2020.nc", var = "temp", orig = "1970-01-01"),
                gdd = list(file = "5_output/1_climate_indices/gdd/fldmean_regions_GDD200_yearly_1763-2020.nc", var = "temp", orig = "1970-01-01"),
                fd = list(file = "5_output/1_climate_indices/fd/fldmean_regions_FD_season_1763-2020.nc", var = "temp", orig = "1970-01-01"),
                wd = list(file ="5_output/1_climate_indices/wd/fldmean_regions_WD_season_1763-2020.nc", var = "precip", orig = "1970-01-01"),
                maesh13d = list(file ="5_output/4_pheno_recon/fldmean_CH_M1_maesh13d_1763-01-01-2020-12-31_2023-06-08_swissplateau.nc", var = "DOY", orig = "1763-01-01"),
                mprua65d = list(file ="5_output/4_pheno_recon/fldmean_CH_M1_mprua65d_1763-01-01-2020-12-31_2023-06-08_swissplateau.nc", var = "DOY", orig = "1763-01-01"),
                mfags13d = list(file ="5_output/4_pheno_recon/fldmean_CH_M1_mfags13d_1763-01-01-2020-12-31_2023-06-07_swissplateau.nc", var = "DOY", orig = "1763-01-01"),
                spring_index = list(file = "5_output/4_pheno_recon/Springindex_phenorecon_2023-06-13.RData"),
                spring_index_mch = list(file = "5_output/4_pheno_recon/Springindex_MCH_2023-06-13.RData"),
                frost_index = list(file = "5_output/4_pheno_recon/CH_frostindex_th0_after_mprua65d_2023-07-07.RData"))

daycounts <- c("cwdi","hwdi","fd","summerdays", "snowdays","wd")
nams <- names(indices)
longnams <- c("Mean temperature °C","Growing degree days","Number of frost day","Number of wet days")


## make a list of the seasonal indices
spring_list <- list()

for(ii in 1:length(indices)){
  print(ii)
  if(!grepl(".RData",indices[[ii]]$file)){
    nc <- nc_open(indices[[ii]]$file)
    dat <- ncvar_get(nc,indices[[ii]]$var)
    time <- as.Date(ncvar_get(nc,"time"),origin = "1970-01-01")
    spring_list[[ii]]  <- list(dat= dat,time = time)
    nc_close(nc)
  } else{
    (load(indices[[ii]]$file)) ##
  }
}

names(spring_list) <- names(indices)[1:7]

spring_list$maesh13d$dat[1] <- NA
spring_list$mprua65d$dat[1] <- NA
spring_list$mfags13d$dat[1] <- NA

### add wet days to 20cr as index
time_wds
index_list_20crv3$index$wd$ind <- array(wds_20cr[iind_wd,jind_wd,], dim = c(4,210))[c(1,3,2,4),]
index_list_20crv3$index$wd$time <- unique(year(time_wds))


#########################################
sigma = 3
lwd = 1.5; lwd2 = 1; lwdabl = 0.5;cx = 1
col = "#115f9a"; col2 = "#76c68f"; col3 = "#c9e52f"; col4 = "#d0f400"

ydates = seq(as.Date("1900-01-01"),as.Date("1900-12-31"), by = "day")
doys = yday(ydates)

abldates <- seq(1760, 2020, 10)
abldatesy <- doys[which(substr(ydates,9,10) %in% c("01"))]

###### spring  #########
png(paste0(plotdir,"/Figure_Timeseries_1763-2020_",Sys.Date(),".png"), width = 2400, height = 300*3, res = 300, pointsize = 6)
par(mfcol = c(3,2), mar = c(0,2,0,2), oma = c(3,2,1,2))

## climate indices
for(ii in 1:length(spring_list[1:4])){
  (svar <- names(spring_list)[[ii]])
  indvals <- spring_list[[svar]]
  if(length(indvals$time) > 258){tind <- month(indvals$time) %in% 4} else {tind = 1:258}
  indvals_sw <- get("index_list")[[svar]]
  indvals_20cr <- index_list_20crv3$index[[svar]]
  indtime <- as.numeric(substr(indvals$time[tind],1,4))
  plot(indtime, indvals$dat[regind,tind],ylab = "", xlab = "", xaxt = "n", col = t_col(col,40), yaxt = "n", type = "l")
  points(indtime, gaussian_wrap(indvals$dat[regind,tind], sigma = sigma), type= "l", lwd = lwd, col = col)
  if(length(indvals_sw) > 0){
    points(indvals_sw[,1], gaussian_wrap(indvals_sw[,ifelse(svar == "gdd",2,"MAM")], sigma = sigma), type= "l", lwd = lwd, col = col2)}
  if(length(indvals_20cr) > 0){
    ifelse(svar== "gdd", indsel <- indvals_20cr$ind[iind,jind,],ifelse(svar == "wd",indsel <- indvals_20cr$ind[3,], indsel <- indvals_20cr$ind[iind,jind,,3]))
    points(indvals_20cr$time, gaussian_wrap(indsel, sigma = sigma), type= "l", lwd = lwd, col = col3)
  }
  mtext(side = 3, line = -1.5, text = paste0(letters[ii],") ",longnams[ii]), adj = 0.01, cex = cx)
  if(svar %in% daycounts){
    axis(side = 2, las = 2)
    #abline(h =  seq(0,40,10), lty = 3, col = "gray75", lwd = lwdabl)
  } else if(svar == "tempmean") {
    axis(side = 2, las = 2, at = seq(1,20,2))
    #abline(h =  seq(1,20,2), lty = 3, col = "gray75", lwd = lwdabl)
  } else{
    axis(side = 2, las = 2, at = seq(0,300,25))
    #abline(h =  seq(1,20,2), lty = 3, col = "gray75", lwd = lwdabl)
  }
  if(ii == 3){axis(side = 1, at = abldates)}
  if(ii == 1){legend("bottom", ncol = 3, col = c(col,col2,col3), legend = c("Swiss recon", "Swiss series", "20crv3"), bty = "n", lwd = lwd, cex = cx)}
}

## pheno indices
days <- seq(as.Date("1900-01-01"),as.Date("1900-07-01"), by = "day")
ind <- substr(days,9,10) %in% c("01","15")
plot(1763:2020, spring_list$maesh13d$dat, type = "l", ylab = "", xlab = "", yaxt = "n", xlim = c(1760, 2020), 
     xaxt = "n",  col = col, lwd = lwd, ylim = c(90,140))
points(1763:2020, spring_list$mfags13d$dat, type = "l",  col = col2, lwd = lwd)
points(1763:2020, spring_list$mprua65d$dat, type = "l",  col = col3, lwd = lwd)
axis(side = 2, at = which(ind) , substr(days,6,10)[ind], las = 2)
#abline(h = which(ind), lty = 2, col = "gray75", lwd = 0.5)
mtext(side = 3, "e) Phenophases", line  = -1.5, cex = cx, adj = 0.01)
legend("bottom", ncol = 4, col = c(col,col2,col3), lwd = lwd, legend = c("Horse chestnut leaf unfolding", "Beech leaf unfolding", "Cherry flowering"), bty = "n", cex =cx)

plot(1763:2020, df_springs_yr$index, type = "l", ylab = "", xlab = "", yaxt = "n", xlim = c(1760, 2020), xaxt = "n",  col = t_col(col,40),ylim = c(-10,24))
points(1763:2020, gaussian_wrap(df_springs_yr$index, sigma = 3), type = "l",  col = col, lwd = lwd)
points(sprg_mch$year, sprg_mch$index, type = "l",  col = t_col(col3, 40))
points(sprg_mch$year, gaussian_wrap(sprg_mch$index, sigma = 3), type = "l",  col = col3, lwd = lwd)
axis(side = 2, las = 2, at = seq(-10,20,10))
#abline(h =  seq(-10,20,5), lty = 2, col = "gray75", lwd = 0.5)
mtext(side = 3, "f) Spring index",line  = -1.5, cex = cx, adj = 0.01)
legend("bottom", ncol = 2, col = c(col,col3), lwd = lwd, legend = c("Reconstruction", "MCH-Index"), bty = "n", cex = cx)
axis(side = 1, at = seq(1760,2020,20) , seq(1760,2020,20), las = 1)

dev.off()


###### spring non smoothed #########
png(paste0(plotdir,"/Figure_Timeseries_1763-2020_v2_",Sys.Date(),".png"), width = 2400, height = 300*3, res = 300, pointsize = 6)
par(mfcol = c(3,2), mar = c(0,2,0,2), oma = c(3,2,1,2))

## climate indices
for(ii in 1:length(spring_list[1:4])){
  (svar <- names(spring_list)[[ii]])
  indvals <- spring_list[[svar]]
  if(length(indvals$time) > 258){tind <- month(indvals$time) %in% 4} else {tind = 1:258}
  indvals_sw <- get("index_list")[[svar]]
  indvals_20cr <- index_list_20crv3$index[[svar]]
  indtime <- as.numeric(substr(indvals$time[tind],1,4))
  plot(indtime, indvals$dat[regind,tind],ylab = "", xlab = "", xaxt = "n", col = t_col(col,40), yaxt = "n", type = "l")
  points(indtime, gaussian_wrap(indvals$dat[regind,tind], sigma = sigma), type= "l", lwd = lwd, col = col)
  if(length(indvals_sw) > 0){
    points(indvals_sw[,1], indvals_sw[,ifelse(svar == "gdd",2,"MAM")], type= "l", lwd = lwd2, col = col2)}
  if(length(indvals_20cr) > 0){
    ifelse(svar== "gdd", indsel <- indvals_20cr$ind[iind,jind,],indsel <- indvals_20cr$ind[iind,jind,,3])
    points(indvals_20cr$time, indsel, type= "l", lwd = lwd2, col = col3)
  }
  points(indtime, indvals$dat[regind,tind],ylab = "", xlab = "", xaxt = "n", col = t_col(col,40), yaxt = "n", type = "l")
  points(indtime, gaussian_wrap(indvals$dat[regind,tind], sigma = sigma), type= "l", lwd = lwd, col = col)
  
  mtext(side = 3, line = -1.5, text = paste0(letters[ii],") ",longnams[ii]), adj = 0.01, cex = cx)
  if(svar %in% daycounts){
    axis(side = 2, las = 2)
  } else if(svar == "tempmean") {
    axis(side = 2, las = 2, at = seq(1,20,2))
  } else{
    axis(side = 2, las = 2, at = seq(0,300,25))
  }
  if(ii == 3){axis(side = 1, at = abldates)}
  if(ii == 1){legend("bottom", ncol = 3, col = c(col,col2,col3), legend = c("Swiss recon", "Swiss series", "20crv3"), bty = "n", lwd = lwd, cex = cx)}
}

## pheno indices
days <- seq(as.Date("1900-01-01"),as.Date("1900-07-01"), by = "day")
ind <- substr(days,9,10) %in% c("01","15")
plot(1763:2020, spring_list$maesh13d$dat, type = "l", ylab = "", xlab = "", yaxt = "n", xlim = c(1760, 2020), 
     xaxt = "n",  col = col, lwd = lwd, ylim = c(90,140))
points(1763:2020, spring_list$mfags13d$dat, type = "l",  col = col2, lwd = lwd)
points(1763:2020, spring_list$mprua65d$dat, type = "l",  col = col3, lwd = lwd)
axis(side = 2, at = which(ind) , substr(days,6,10)[ind], las = 2)
#abline(h = which(ind), lty = 2, col = "gray75", lwd = 0.5)
mtext(side = 3, "e) Phenophases", line  = -1.5, cex = cx, adj = 0.01)
legend("bottom", ncol = 4, col = c(col,col2,col3), lwd = lwd, legend = c("Horse chestnut leaf unfolding", "Beech leaf unfolding", "Cherry flowering"), bty = "n", cex =cx)

## spring index
plot(1763:2020, df_springs_yr$index, type = "l", ylab = "", xlab = "", yaxt = "n", xlim = c(1760, 2020), xaxt = "n",  col = t_col(col,40),ylim = c(-10,24))
points(1763:2020, gaussian_wrap(df_springs_yr$index, sigma = 3), type = "l",  col = col, lwd = lwd)
points(sprg_mch$year, sprg_mch$index, type = "l",  col = t_col(col3, 40))
points(sprg_mch$year, gaussian_wrap(sprg_mch$index, sigma = 3), type = "l",  col = col3, lwd = lwd)
axis(side = 2, las = 2, at = seq(-10,20,10))
#abline(h =  seq(-10,20,5), lty = 2, col = "gray75", lwd = 0.5)
mtext(side = 3, "f) Spring index",line  = -1.5, cex = cx, adj = 0.01)
legend("bottom", ncol = 2, col = c(col,col3), lwd = lwd, legend = c("Reconstruction", "MCH-Index"), bty = "n", cex = cx)
axis(side = 1, at = seq(1760,2020,20) , seq(1760,2020,20), las = 1)

dev.off()

