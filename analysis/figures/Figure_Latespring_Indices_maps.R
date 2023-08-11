##### PLOT LATE SPRING INDICES #######
rm(list=ls())

### packages
library(ncdf4);library(lubridate);library(RColorBrewer)

### helpfer functions
source("R/image_funs.R")

### directores
plotdir <- "manuscript/figures/"

### topography
altth <- 1600
nc <- nc_open("../swiss_arm/code/data_input/topo.swiss02_ch01r.swisscors.nc")
topo <- ncvar_get(nc, "height")
topo_help <- topo
topo_help[topo >= altth] <- 10
topo_help[topo < altth] <- 1
topo_N <-ncvar_get(nc, "chy") + 1000000
topo_E <-ncvar_get(nc, "chx") + 2000000
nc_close(nc)

springs <- c(1770,1785,1816,1837,1853,1932)

### indices
indices <- list(fd = list(file = "data/01_climate_indices/fd/FD_season_1763-2020.nc", var = "temp", orig = "1970-01-01"),
                mprua65d = list(file ="data/04_pheno_recon/CH_PTTs_mprua65d_1763-01-01-2020-12-31_2023-08-02.nc", var = "DOY", orig = "1763-01-01"),
                snow = list(file ="data/01_climate_indices/snowdays/SNOWDAYS_season_1763-2020.nc", var = "temp", orig = "1970-01-01"),
                wetdays = list(file ="data/01_climate_indices/wd/WD_season_1763-2020.nc", var = "precip", orig = "1970-01-01"),
                tempmean = list(file = "data/01_climate_indices/tempmean/TEMPMEAN_season_1763-2020.nc", var = "temp", orig = "1970-01-01"))

for(ii in 1:length(indices)){
  (var <- names(indices)[ii])
  nc <-nc_open(indices[[ii]]$file)
  time <- as.Date(ncvar_get(nc,"time"), origin = indices[[ii]]$orig)
  tind <- which(year(time) %in% springs )
  if(length(tind) > 10){tind <- which(year(time) %in% springs  & month(time) %in% c(3:5))}
  dat <- ncvar_get(nc,indices[[ii]]$var)[,,tind]
  dat[dat < -888] <- NA
  assign(paste0(var,"_mat"),dat)
  assign(paste0(var,"_time"),time[tind])
  nc_close(nc)
}

### calculate climatology
for(ii in 1:length(indices)){
  (var <- names(indices)[ii])
  nc <-nc_open(indices[[ii]]$file)
  time <- as.Date(ncvar_get(nc,"time"), origin = indices[[ii]]$orig)
  if(var %in% c("fd","snow","wetdays","tempmean")){
    tind <- which(year(time) %in% 1871:1900  & month(time) %in% c(3:5))
    assign(paste0(var,"_mat_clim"),apply(ncvar_get(nc,indices[[ii]]$var)[,,tind],1:2,mean))
  }
  nc_close(nc)
}

### mean pheno clim
for(ii in 2:5){
  (var <- names(indices)[ii])
  nc <-nc_open(indices[[ii]]$file)
  time <- as.Date(ncvar_get(nc,"time"), origin = indices[[ii]]$orig)
  tind <- which(year(time) %in% 1871:1900)
  assign(paste0(var,"_help"),ncvar_get(nc,indices[[ii]]$var)[,,tind])
  nc_close(nc)
}

mprua65d_clim <- apply(mprua65d_help ,1:2,mean)

nc <-nc_open(indices[[1]]$file)
lat <- ncvar_get(nc,"N")
lon <- ncvar_get(nc,"E")
nc_close(nc)


############## the three latest springs ####################
lty = 1; adj = 0.01
cx = 1; marscale = c(3,2,1,2); mline = 2.5
lsprgs <- c(1785,1837,1853)

png(paste0(plotdir,"Figure_Latesprings_Indices_",Sys.Date(),".png"), width = 3600, height = 2000, res = 300, pointsize = 7)
layout(matrix(1:((length(lsprgs) + 1)*5), nrow = 4, ncol = 5, byrow = F), widths = c(1), heights = c(1,1,1,0.22))
par(oma=c(1,1,1,2))
### Temp
par(mar=c(1,1,1,1))
ind <-  "tempmean"
brksp <- seq(-3.5,3.5,0.5)
colp <-  rev(colorRampPalette(brewer.pal(n =5, name = "RdBu"))(length(brksp)-1))
for(l in 1:length(lsprgs)){
  mat <- get(paste0(ind,"_mat"))[,,which(springs == lsprgs[l])] - get(paste0(ind,"_mat_clim"))
  mat[mat < brksp[1]] <- brksp[1]
  mat[mat > brksp[length(brksp)]] <- brksp[length(brksp)]
  image(lon,lat,mat, bty = "n", yaxt = "n", xaxt = "n", col = colp, breaks = brksp)
  mtext(side = 2, text = lsprgs[l])
  if(l == 1) mtext(side = 3, text = "Mean temperature")
}
image_scale(breaks=brksp,col=colp, axis.pos=1, scale_lab="°C",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=1,ats=brksp[seq(1,length(brksp),2)],labs=brksp[seq(1,length(brksp),2)], cex=cx,axis.pos2=1,line=mline, adj=0))


### Pheno
par(mar=c(1,1,1,1))
brksp <- seq(-25,25,5)
colp <- colorRampPalette(c("#0E6134","#2BAC66","#ABCF63","#E8F69E","#D7E3EE","#8FB3FF","#7F97FF","#0570B0","#023858"))(length(brksp)-1)
for(l in 1:length(lsprgs)){
  mat <- mprua65d_mat[,,which(springs == lsprgs[l])] - mprua65d_clim
  mat[mat < brksp[1]] <- brksp[1]
  mat[mat > brksp[length(brksp)]] <- brksp[length(brksp)]
  mat[topo > altth] <- NA
  image(lon,lat,topo_help, bty = "n", yaxt = "n", xaxt = "n", col = c("gray90","gray65"), breaks = c(0,1,10), xlab = "", ylab = "")
  image(lon,lat,mat, bty = "n", yaxt = "n", xaxt = "n", col = colp, breaks = brksp, xlab = "", ylab = "", add  = T)
  #contour(lon,lat,meanpheno_mat_clim, bty = "n", yaxt = "n", xaxt = "n", levels = c(120), add = T, col = "gray30", drawlabels = FALSE, lwd = 0.5)
  #contour(lon,lat,meanpheno_mat_clim, bty = "n", yaxt = "n", xaxt = "n", levels = c(140), add = T, col = "gray40", drawlabels = FALSE, lwd = 0.5)
  if(l == 1) mtext(side = 3, text = "Cherry flowering")
}
image_scale(breaks=brksp,col=colp, axis.pos=1, scale_lab="Doy",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=1,ats=brksp[seq(1,length(brksp),2)],labs=brksp[seq(1,length(brksp),2)], cex=cx,axis.pos2=1,line=mline, adj=0))

### Frost
par(mar=c(1,1,1,1))
ind <- "fd"
brksp <- seq(-30,30,5)
colp <- colorRampPalette(c("#0E6134","#2BAC66","#ABCF63","#E8F69E","#D7E3EE","#B5CAFF","#8FB3FF","#7F97FF","#0570B0","#023858"))(length(brksp)-1)
for(l in 1:length(lsprgs)){
  mat <- get(paste0(ind,"_mat"))[,,which(springs == lsprgs[l])] - get(paste0(ind,"_mat_clim"))
  mat[mat < brksp[1]] <- brksp[1]
  mat[mat > brksp[length(brksp)]] <- brksp[length(brksp)]
  image(lon,lat,mat, bty = "n", yaxt = "n", xaxt = "n", col = colp, breaks = brksp)
  if(l == 1) mtext(side = 3, text = "Frost days")
}
image_scale(breaks=brksp,col=colp, axis.pos=1, scale_lab="Days",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=1,ats=brksp[seq(1,length(brksp),2)],labs=brksp[seq(1,length(brksp),2)], cex=cx,axis.pos2=1,line=mline, adj=0))

### Snow
par(mar=c(1,1,1,1));ind <- "snow"
brksp <- seq(-25,25,5)
colp <- colorRampPalette(c("#0E6134","#2BAC66","#ABCF63","#E8F69E","#D7E3EE","#8FB3FF","#7F97FF","#0570B0","#023858"))(length(brksp)-1)
for(l in 1:length(lsprgs)){
  mat <- get(paste0(ind,"_mat"))[,,which(springs == lsprgs[l])] - get(paste0(ind,"_mat_clim"))
  mat[mat < brksp[1]] <- brksp[1]
  mat[mat > brksp[length(brksp)]] <- brksp[length(brksp)]
  image(lon,lat,mat, bty = "n", yaxt = "n", xaxt = "n", col = colp, breaks = brksp)
  if(l == 1) mtext(side = 3, text = "Snow days")
}
image_scale(breaks=brksp,col=colp, axis.pos=1, scale_lab="Days",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=1,ats=brksp[seq(1,length(brksp),2)],labs=brksp[seq(1,length(brksp),2)], cex=cx,axis.pos2=1,line=mline, adj=0))

### Wet
par(mar=c(1,1,1,1))
ind <- "wetdays"
brksp <- seq(-25,25,5)
colp <- colorRampPalette(c("peru","gray98","cornflowerblue"))(length(brksp)-1)
for(l in 1:length(lsprgs)){
  mat <- get(paste0(ind,"_mat"))[,,which(springs == lsprgs[l])] - get(paste0(ind,"_mat_clim"))
  mat[mat < brksp[1]] <- brksp[1]
  mat[mat > brksp[length(brksp)]] <- brksp[length(brksp)]
  image(1:nrow(mat),1:ncol(mat),mat, bty = "n", yaxt = "n", xaxt = "n", col = colp, breaks = brksp)
  if(l == 1) mtext(side = 3, text = "Wet days")
}
image_scale(breaks=brksp,col=colp, axis.pos=1, scale_lab="Days",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=1,ats=brksp[seq(1,length(brksp),2)],labs=brksp[seq(1,length(brksp),2)], cex=cx,axis.pos2=1,line=mline, adj=0))

dev.off()

