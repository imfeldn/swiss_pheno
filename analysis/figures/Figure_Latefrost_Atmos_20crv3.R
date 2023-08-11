###### LATE FROST MAPS WITH 20CRV3 ######
rm(list=ls())

### packages
library(ncdf4);library(lubridate);library(RColorBrewer)

### helpfer functions
source("R/helpfuns.R")
source("R/image_funs.R")

### directories
plotdir <- "manuscript/figures/"

## read files
t2m_file <- "../swiss_indices/2_data/20crv3/t2m/t2m_daily_anom_1806-2015.nc"
gph_file <- "../swiss_indices/2_data/20crv3/gph500/gph_daily_1806-2015.nc"

ncg <- nc_open(gph_file)
gtime <- as.Date(ncvar_get(ncg,"time")/24, origin = "1800-01-01")
glon <- ncvar_get(ncg,"longitude")
glon[glon > glon[length(glon)]] <- glon[glon > glon[length(glon)]] -360
glat <- ncvar_get(ncg,"latitude");glat <- rev(glat)
yind1 <- which(year(gtime) == 1873);yind2 <- which(year(gtime) == 1957)
gph1 <- ncvar_get(ncg,"gh", c(1,1,1,yind1[1]), c(-1,-1,1,length(yind1)))
gph2 <- ncvar_get(ncg,"gh", c(1,1,1,yind2[1]), c(-1,-1,1,length(yind2)))
gph <- abind::abind(gph1,gph2)[,length(glat):1,];rm(gph1,gph2)
gtimesel <- gtime[c(yind1, yind2)]
nc_close(ncg)

nct <- nc_open(t2m_file)
ttime <- as.Date(ncvar_get(nct,"time")/24, origin = "1800-01-01")
tlon <- ncvar_get(nct,"longitude")
tlon[tlon > tlon[length(tlon)]] <- tlon[tlon > tlon[length(tlon)]] -360
tlat <- ncvar_get(nct,"latitude");tlat <- rev(tlat)
yind1 <- which(year(ttime) == 1873);yind2 <- which(year(ttime) == 1957)
t2m1 <- ncvar_get(nct,"t2m", c(1,1,yind1[1]), c(-1,-1,length(yind1)))
t2m2 <- ncvar_get(nct,"t2m", c(1,1,yind2[1]), c(-1,-1,length(yind2)))
t2m <- abind::abind(t2m1,t2m2)[,length(tlat):1,];rm(t2m1,t2m2)
timesel <- ttime[c(yind1, yind2)]
nc_close(nct)

brks <- seq(-10,10,2)#seq(100,145,2)
cols <- rev(colorRampPalette(brewer.pal(n = 5, name = "RdBu"))(length(brks)-1))

### plot March accumulation and frost day
marscale <- c(1,0.5,1.5,2); cx = 1
png(paste0(plotdir,"Figure_Latefrosts_Atmos.png"), width = 2900, height = 800, res = 300, pointsize = 8)
layout(matrix(1:5, nrow = 1, ncol = 5, byrow = T), heights = c(1,1), widths = c(rep(1,4),0.18))
par(mar = c(1,0.5,1.5,0.5), oma = c(2,2,1,1.5))

for(ii in c(1873,1957)){
  tind <- year(timesel) == ii & month(timesel) %in% 3
  image(tlon,tlat,apply(t2m[,,tind],1:2,mean), xlim = c(-20,50), breaks = brks, col = cols, xaxt = "n", yaxt = "n")
  maps::map("world",add = T, col = "gray55")
  maps::map("world",region = "Switzerland", add = T, col = "red")
  contour(glon,glat, apply(gph[,,tind],1:2,mean), add = T, col = "gray45")
  #if(ii == 1873) mtext(side = 3, adj = 0.01, line = 0.1, "b)")
  mtext(side = 3, line = 0.1, text = paste0("March ", ii))
}

for(ii in 1:2){
  tind <- as.character(timesel) == as.character(c("1873-04-26","1957-05-08")[ii])
  dat <- t2m[,,tind]
  dat[dat < brks[1]] <- brks[1]
  dat[dat > brks[length(brks)]] <- brks[length(brks)]
  image(tlon,tlat,dat, xlim = c(-20,50), breaks = brks, col = cols, xaxt = "n", yaxt = "n")
  maps::map("world",add = T, col = "gray55")
  maps::map("world",region = "Switzerland", add = T, col = "red")
  contour(glon,glat, apply(gph[,,tind],1:2,mean), add = T, col = "gray45")
  mtext(side = 3, line = 0.1, text = as.character(c("1873-04-26","1957-05-08")[ii]))
}

image_scale(breaks=brks,col=cols, axis.pos=4, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=4,ats=brks[seq(1,length(brks),2)],labs=brks[seq(1,length(brks),2)], cex=cx,axis.pos2=4,line=3, adj=-0.07))
dev.off()

