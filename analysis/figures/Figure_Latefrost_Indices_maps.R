########  PLOT LATEFROST INDICES #####
rm(list=ls())

### packages
library(ncdf4);library(lubridate)
library(viridis);library(RColorBrewer)

### functions
source("R/helpfuns.R")
source("R/image_funs.R")

### directories
plotdir <- "manuscript/figures/"

### load altitude
altth <- 1600
nc <- nc_open("../swiss_arm/code/data_input/topo.swiss02_ch01r.swisscors.nc")
topo <- ncvar_get(nc, "height")
topobl <- which(topo < altth)
topo_help <- topo
topo_help[topo >= altth] <- 10
topo_help[topo < altth] <- 1
topo_N <-ncvar_get(nc, "chy") + 1000000
topo_E <-ncvar_get(nc, "chx") + 2000000
nc_close(nc)

### get wgs coordinates
nchelp <- nc_open("../swiss_arm/code/data_input/lonlats_temp.nc")
lon <- ncvar_get(nchelp,"lon")[,1]
lat <- ncvar_get(nchelp,"lat")[1,]

pv <- "mprua65d"
indices <- list(tmean_after_pheno = list(file = paste0("data/04_pheno_recon/CH_frostindex_th0_after_mprua65d_5daysprior_2023-08-03.RData")),
                pheno = list(file =paste0("data/04_pheno_recon/CH_PTTs_mprua65d_1763-01-01-2020-12-31_2023-08-02.nc"), var = "DOY", orig = "1763-01-01"),
                lastfrost = list(file ="data/01_climate_indices/lfd/LFD_yearly_1763-2020.nc", var = "temp", orig = "1970-01-01"))

for(ii in 1:length(indices)){
  if(grepl("RData",indices[[ii]]$file)){
    (var <- names(indices)[ii])
    (load(indices[[ii]]$file))
    assign(paste0(var,"_mat"),tmean_after_pheno)
  } else{
    (var <- names(indices)[ii])
    nc <-nc_open(indices[[ii]]$file)
    assign(paste0(var,"_time"),as.Date(ncvar_get(nc,"time"), origin = indices[[ii]]$orig))
    assign(paste0(var,"_mat"),ncvar_get(nc,indices[[ii]]$var))
  }
}

tmean_after_pheno_mat[,,1] <- NA

## calculate climatology
for(ii in 2:length(indices)){
   (var <- names(indices)[ii])
   nc <-nc_open(indices[[ii]]$file)
   time <- as.Date(ncvar_get(nc,"time"), origin = indices[[ii]]$orig)
   tind <- which(year(time) %in% 1871:1900)
   assign(paste0(var,"_mat_clim"),apply(ncvar_get(nc,indices[[ii]]$var)[,,tind],1:2,mean))
}

### calc time series for area below 2000masl?
topo_mask <- (topo < altth)
frostindex_ts <- sapply(1:dim(tmean_after_pheno_mat)[3], function(x){mean(tmean_after_pheno_mat[,,x][topo_mask], na.rm = T)})
frostarea  <- sapply(1:dim(tmean_after_pheno_mat)[3], function(x){sum(!is.na(tmean_after_pheno_mat[,,x][topo_mask]), na.rm = T)})

df_frost <- data.frame(year =1763:2020, frostindex_ts, frostarea)
View(df_frost)

# png(paste0(plotdir,"frostindex_area_mean_bl",altth,"_",Sys.Date(),".png"), width = 1600)
# par(mfrow= c(2,1), mar = c(2,2,1,1))
# plot(1763:2020, abs(df_frost$frostindex_ts), type = "l", lwd = 2);abline(v = seq(1760,2020,5), lty = 3, col = "gray75");title("frostindex", line = -1.5)
# plot(1763:2020, df_frost$frostarea, type = "l", lwd = 2);abline(v = seq(1760,2020,5), lty = 3, col = "gray75");title("affected area bl 1600", line = -1.5)
# dev.off()

lty <- 1; adj <- 0.01
cx <- 1.3; marscale <- c(1,3,2,3)
selyrs <- c(1873,1957)

days <- substr(seq(as.Date("1900-01-01"),as.Date("1900-07-01"), by = "day"),6,10)
spind <- 10

png(paste0(plotdir,"Figure_Latefrosts_Indices_",Sys.Date(),".png"), width = 2300, height = 1300, res = 300, pointsize = 10)
layout(matrix(1:9, nrow = 3, ncol = 3, byrow = F), widths = 1, heights = c(1,1,0.22))
par(oma = c(2,2,2,1))
par(mar=c(1,1,1,1));ind <- names(indices)[1]
brksp <-seq(-14,0,2)
colp <- rev(colorRampPalette(c("lightskyblue","darkblue"))(length(brksp)-1))
for(l in 1:2){
  mat <- get(paste0(ind,"_mat"))[,,which(1763:2020 == selyrs[l])]
  mat[mat < brksp[1]] <- brksp[1]
  mat[mat > brksp[length(brksp)]] <- brksp[length(brksp)]
  mat[topo > altth] <- NA
  image(lon,lat,topo_help, bty = "n", yaxt = "n", xaxt = "n", col = c("gray90","gray45"), breaks = c(0,1,10), xlab = "", ylab = "")
  image(lon,lat,mat, bty = "n", yaxt = "n", xaxt = "n", col = colp, breaks = brksp, add = T)
  mtext(side = 2, text = selyrs[l])
  if(l == 1) mtext(side = 3, text = "Frost index")
}
image_scale(breaks=brksp,col=colp, axis.pos=1, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=1,ats=brksp[seq(1,length(brksp),2)],labs=brksp[seq(1,length(brksp),2)], cex=cx,axis.pos2=1,line=3, adj=-0.07))

par(mar=c(1,1,1,1))
brksp <- seq(85,150,5)
colp <- colorRampPalette(c("#0E6134","#2BAC66","#ABCF63","#E8F69E","#D7E3EE","#8FB3FF","#7F97FF","#0570B0","#023858"))(length(brksp)-1)
for(l in 1:2){
  mat <- lastfrost_mat[,,which(1763:2020 == selyrs[l])]
  mat[mat < brksp[1]] <- brksp[1]
  mat[mat > brksp[length(brksp)]] <- brksp[length(brksp)]
  mat[topo > 1600] <- NA
  image(lon,lat,topo_help, bty = "n", yaxt = "n", xaxt = "n", col = c("gray90","gray45"), breaks = c(0,1,10), xlab = "", ylab = "")
  image(lon,lat,mat, col = colp, breaks = brksp, add = T)
  anom <- mat - lastfrost_mat_clim
  spaceind <- which(anom > 15, arr.ind = T)
  points(lon[spaceind[which(spaceind[,1] %in%  seq(1,length(lon),10)),1]],lat[spaceind[which(spaceind[,1] %in%  seq(1,length(lon),10)),2]],
         pch = 20, cex = 0.0001, col = t_col("magenta", 100))

  if(l == 1) mtext(side = 3, text = "Last frost day")
}
labs <- days[brksp[seq(1,length(brksp),2)]]

image_scale(breaks=brksp,col=colp, axis.pos=1, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=1,ats=brksp[seq(1,length(brksp),2)],labs=labs,
                            cex=cx,axis.pos2=1,line=3, adj=-0.07))

par(mar=c(1,1,1,1));ind <- names(indices)[2]
for(l in 1:2){
  mat <- get(paste0(ind,"_mat"))[,,1763:2020 == selyrs[l]]
  mat[mat < brksp[1]] <- brksp[1]
  mat[mat > brksp[length(brksp)]] <- brksp[length(brksp)]
  mat[topo > 1600] <- NA
  image(lon,lat,topo_help, bty = "n", yaxt = "n", xaxt = "n", col = c("gray90","gray45"), breaks = c(0,1,10), xlab = "", ylab = "")
  image(lon,lat,mat , bty = "n", col = colp, breaks = brksp, add = T)
  anom <- mat - pheno_mat_clim
  spaceind <- which(anom < -15, arr.ind = T)
  points(lon[spaceind[which(spaceind[,1] %in%  seq(1,length(lon),10)),1]],lat[spaceind[which(spaceind[,1] %in%  seq(1,length(lon),10)),2]],
         pch = 20, cex = 0.0001, col = t_col("magenta", 100))
  if(l == 1) mtext(side = 3, text = "Cherry flowering")
}
image_scale(breaks=brksp,col=colp, axis.pos=1, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=1,ats=brksp[seq(1,length(brksp),2)],labs=labs, cex=cx,axis.pos2=1,line=3, adj=-0.07))

dev.off()
