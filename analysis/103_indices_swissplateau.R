#### CALCULATE INDICES FOR THE SWISS PLATEAU SERIES FROM BRUGNARA ET AL. 2022 #######
### https://doi.org/10.5194/cp-18-2357-2022

### packages
library(ncdf4);library(lubridate);library(mmand)

### help function
source("R/helpfuns.R")
source("R/helpfuns_swissplateau.R")
source("R/helpfuns_temp_quantiles.R")

### define indices
indices <- c("gsstart","gsend","gsl","ffd","lfd","fd","gdd","summerdays","wsdi","csdi","tempmean")

### read data
merge <- F
if(merge){
  mch_dat <- read.table("../swiss_arm/data/mch/nbcn/temp/order_99106_data.txt", header = T, skip = 2)
  ber_dat <- mch_dat[mch_dat$stn == "BER",]
  sma_dat <- mch_dat[mch_dat$stn == "SMA",]
  mix <-  data.frame(Date = as.Date(ber_dat$time,format = "%Y%m%d"),Ta_mean = (as.numeric(ber_dat$ths200d0) + as.numeric(sma_dat$ths200d0))/2)
  sw_dat <- read.table("../swiss_arm/data/EarlyInstrumentalTemperature/SwissPlateau/SwissPlateau_1756-1863_daily.csv", header = T, sep = ",")
  mondiff <- c(-0.8500,-0.7750,-0.8125,-0.8875,-1.1125,-1.0250,-0.9250,-0.8500,-0.7750,-0.7625,-0.8375,-0.7750) ## these are the monthly homogenization values define by Y.Brugnara
  sw_dat_adj <- sw_dat$Ta_mean + mondiff[as.numeric(substr(sw_dat$Date,6,7))]
  sw_dat <- data.frame(Date = as.Date(sw_dat$Date),Ta_mean = sw_dat_adj)
  sw_full <- rbind(sw_dat, mix)
  save(sw_full, file = "data/05_swissplateau/swissplateau_daily_hom_1756-2021.RData")
} else{
  (load("data/05_swissplateau/swissplateau_daily_hom_1756-2021.RData"))
}

index_list <- list()
ii <- 11

for(ii in 1:length(indices)){
  print(indices[ii])
  func <- get(indices[ii])
  if(indices[ii] == "gdd"){
    index_list[[ii]] <- func(dat = sw_full$Ta_mean, date = sw_full$Date, th = 200)
  } else if(indices[ii] %in% c("csdi","wsdi")){
    qth <- ifelse(indices[ii] == "cwdi", 0.2, 0.8)
    qdat <- calc_quants_station(temp = sw_full$Ta_mean, date = sw_full$Date, qth = qth)
    index_list[[ii]] <- func(temp = sw_full$Ta_mean, dates = sw_full$Date, qdat = qdat)
  } else {
    index_list[[ii]] <- func(dat = sw_full$Ta_mean, date = sw_full$Date)
  }
}


### add precipitation indices as comparison
### calculate wet days and snow days individually and then create average
mch_dat <- read.table("../swiss_arm/data/mch/nbcn/temp/order_99106_data.txt", header = T, skip = 2)
ber_temp <- mch_dat[mch_dat$stn == "BER",]
sma_temp <- mch_dat[mch_dat$stn == "SMA",]
sma_temp_dates <- as.Date(paste(substr(mch_dat[mch_dat$stn == "SMA",2],1,4),
                    substr(mch_dat[mch_dat$stn == "SMA",2],5,6),
                    substr(mch_dat[mch_dat$stn == "SMA",2],7,8), sep = "-"))
ber_temp_dates <- as.Date(paste(substr(mch_dat[mch_dat$stn == "BER",2],1,4),
                        substr(mch_dat[mch_dat$stn == "BER",2],5,6),
                        substr(mch_dat[mch_dat$stn == "BER",2],7,8), sep = "-"))

mch_dat <- read.table("../swiss_arm/data/mch/nbcn/precip/order_99111_data.txt", header = T, skip = 2)
ber_precip <- mch_dat[mch_dat$stn == "BER",]
sma_precip <- mch_dat[mch_dat$stn == "SMA",]
sma_precip_dates <- paste(substr(mch_dat[mch_dat$stn == "SMA",2],1,4),
                        substr(mch_dat[mch_dat$stn == "SMA",2],5,6),
                        substr(mch_dat[mch_dat$stn == "SMA",2],7,8), sep = "-")
ber_precip_dates <- paste(substr(mch_dat[mch_dat$stn == "BER",2],1,4),
                        substr(mch_dat[mch_dat$stn == "BER",2],5,6),
                        substr(mch_dat[mch_dat$stn == "BER",2],7,8), sep = "-")

### check if dates are all equal
all(ber_precip_dates == ber_temp_dates)
all(sma_precip_dates == sma_temp_dates)
all(sma_precip_dates == ber_temp_dates)

## calc wet days
seas_facts <- get_seas_facts(ber_temp_dates)

ber_wd <- (as.numeric(ber_precip$rhs150d0) >= 1) + 0
sma_wd <- (as.numeric(sma_precip$rhs150d0) >= 1) + 0

ber_tg2 <- (as.numeric(ber_temp$ths200d0) < 2) + 0
sma_tg2 <- (as.numeric(sma_temp$ths200d0) < 2) + 0

snow_ber <- ((ber_wd + ber_tg2) == 2) + 0
snow_sma <- ((sma_wd + sma_tg2) == 2) + 0

ber_snow_seas <- aggregate(snow_ber, list(seas_facts), sum_Xpercna, na = 0.1)
sma_snow_seas <- aggregate(snow_sma, list(seas_facts), sum_Xpercna, na = 0.1)
ber_wd_seas <- aggregate(ber_wd, list(seas_facts), sum_Xpercna, na = 0.1)
sma_wd_seas <- aggregate(sma_wd, list(seas_facts), sum_Xpercna, na = 0.1)
yrs <- as.numeric(substr(ber_snow_seas$Group.1,1,4))

ber_snow_seas <- ber_snow_seas[yrs < 2021,];sma_snow_seas <- sma_snow_seas[yrs < 2021,]
ber_wd_seas <- ber_wd_seas[yrs < 2021,];sma_wd_seas <- sma_wd_seas[yrs < 2021,]

sw_snow <- (ber_snow_seas$x + sma_snow_seas$x)/2
sw_wd <- (ber_wd_seas$x + sma_wd_seas$x)/2

uyrs <- unique(yrs[yrs < 2021])
sw_snow_mat <-cbind(uyrs, matrix(sw_snow, ncol = 4, nrow = length(uyrs), byrow = T))
colnames(sw_snow_mat) <- c("year", "DJF", "JJA","MAM","SON")

sw_wd_mat <-cbind(uyrs, matrix(sw_wd, ncol = 4, nrow = length(uyrs), byrow = T))
colnames(sw_wd_mat) <- c("year", "DJF", "JJA","MAM","SON")

index_list[[ii + 1]] <- as.data.frame(sw_snow_mat)
index_list[[ii + 2]] <- as.data.frame(sw_wd_mat)

names(index_list) <- c(indices,"snowdays","wd")
save(index_list, file = paste0("data/01_climate_indices/swissplateau_stations_indices_",Sys.Date(),".RData"))

