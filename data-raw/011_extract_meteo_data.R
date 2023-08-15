####### EXTRACT TEMPERATURE DATA FROM CLOSEST GRID CELLS FOR SELECTED SERIES #############

rm(list=ls())

#### packages
library(ncdf4)
library(stats)
library(geosphere)

#### directories
indir_nc <- "../swiss_recon/temp/"

#### select classes
classes <- "123_longseries"

yrs <- 1951:2020
nyr <- length(yrs)

### read grid data
ft <- nc_open(paste(indir_nc,"CH_temp_EnKF_1951-01-01-1951-12-31.nc",sep="")) # nonly for coordinates
lon <- ft$dim[[2]]$vals
lat <- ft$dim[[3]]$vals
nc_close(ft)

### Umrechnen Gitterdaten in WGS84
londegb <- matrix(NA,nrow=length(lon),ncol=length(lat))
latdegb <- matrix(NA,nrow=length(lon),ncol=length(lat))
Ep <- (lon-2600000)/1000000
Np <- (lat-1200000)/1000000

for (i in 1:length(lon)){for(j in 1:length(lat)){
  londegb[i,j] <- (2.6779094+4.728982*Ep[i]+0.791484*Ep[i]*Np[j]+0.1306*Ep[i]*(Np[j]^2)-0.0436*(Ep[i]^3))*100/36
  latdegb[i,j] <- (16.9023892+3.238272*Np[j]-0.270978*(Ep[i]^2)-0.002528*(Np[j]^2)-0.0447*(Ep[i]^2)*Np[j]-0.0140*(Np[j]^3))*100/36
}}


### list phenophases for extraction
phenovals <- c( "mprua65d", "maesh13d")#"mfags13d","mcora65d","mcora13d","manen65d","mtaro65d","mtusf65d","mcarp65d","mlard13d")
pv <- phenovals[1]

for(pv in phenovals){

  print(pv)

  load(paste0("data/02_pheno_net/pheno_meta_",pv,"_class",classes,".RData"))
  load(paste0("data/02_pheno_net/stations_",pv,"_class",classes,".RData"))

  nstat <- nrow(station_list_sel)
  obslon <- station_list_sel$Laengengrad
  obslat <- station_list_sel$Breitengrad
  obsnam <- as.character(station_list_sel$Abk.)
  obsalt <- station_list_sel$Stationshoehe

  ################################################
  ### get temperature series for every station
  ################################################
  sellon <- rep(NA,nstat)
  sellat <- rep(NA,nstat)

  for (a in 1:nstat){
    dist <- (((obslon[a]-londegb)^2)+((0.682*(obslat[a]-latdegb))^2)^0.5)
    sellon[a] <- which(dist==min(dist),arr.ind=T)[1]
    sellat[a] <- which(dist==min(dist),arr.ind=T)[2]
  }

  ### extract station temperatures
  meteo_dat <- array(NA,dim=c(nstat,nyr + 1,365))
  startyr <- 1950
  tind_yrend2 <- 264:365

  for (yr in startyr:2020){
    print(yr)

    if (yr>1960){

      ft <- nc_open(paste(indir_nc,"CH_temp_TabsD_",yr,".nc",sep=""))
      ty <- ncvar_get(ft,varid="TabsD")
      tind_yrend <- (dim(ty)[3] - 101):dim(ty)[3]
      tind_start <- 1:263

    for (a in 1:nstat){
      meteo_dat[a,yr - startyr + 1,(length(tind_yrend)+1):365] <- round(ty[sellon[a],sellat[a],tind_start],2) # write the data from the current spring to doy 102 from current year
      if(yr < 2020) meteo_dat[a,yr - startyr + 2,1:length(tind_yrend)] <- round(ty[sellon[a],sellat[a],tind_yrend],2) ## write the data from the current winter to doy 1 of next year
    }
    nc_close(ft)
    } else {

      ft <- nc_open(paste(indir_nc,"CH_temp_EnKF_",yr,"-01-01-",yr,"-12-31.nc",sep=""))
      ty <- ncvar_get(ft,varid="temp")
      tind_yrend <- (dim(ty)[3] - 101):dim(ty)[3]
      tind_start <- 1:263

      for (a in 1:nstat){
        meteo_dat[a,yr - startyr + 1,(length(tind_yrend)+1):365] <- round(ty[sellon[a],sellat[a],tind_start],2) # write the data from the current spring to doy 102 from current year
        meteo_dat[a,yr - startyr + 2,1:length(tind_yrend)] <- round(ty[sellon[a],sellat[a],tind_yrend],2) ## write the data from the current winter to doy 1 of next year
      }

      nc_close(ft)
    }
  }

  # remove first year since it is not needed
  meteo_dat <- meteo_dat[,2:dim(meteo_dat)[2],]
  rownames(meteo_dat) <- obsnam

  saveRDS(
    meteo_dat,
    file = paste0("data/02_pheno_net/meteo_",pv,"_class",classes,".rds"),
    compress = "xz"
    )
}
