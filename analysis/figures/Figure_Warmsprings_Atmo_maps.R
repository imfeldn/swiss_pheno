###### late frost maps with 20crv3 ######
rm(list=ls())

### packages
library(ncdf4);library(lubridate);library(RColorBrewer)

### functions
source("R/helpfuns.R")
source("R/image_funs.R")

### directories
plotdir <- "manuscript/figures/"

## read files
t2m_file <- "../swiss_indices/2_data/20crv3/t2m/t2m_daily_anom_1806-2015.nc"
gph_file <- "../swiss_indices/2_data/20crv3/gph500/gph_daily_1806-2015.nc"
gph_anom_file <- "../swiss_indices/2_data/20crv3/gph500/gph_daily_anom_1806-2015.nc"
blocking_files <- ""

ncg <- nc_open(gph_file)
gtime <- as.Date(ncvar_get(ncg,"time")/24, origin = "1800-01-01")
glon <- ncvar_get(ncg,"longitude")
glon[glon > glon[length(glon)]] <- glon[glon > glon[length(glon)]] -360
glat <- ncvar_get(ncg,"latitude");glat <- rev(glat)
yind1 <- which(year(gtime) == 1862);yind2 <- which(year(gtime) == 2011)
gph1 <- ncvar_get(ncg,"gh", c(1,1,1,yind1[1]), c(-1,-1,1,length(yind1)))
gph2 <- ncvar_get(ncg,"gh", c(1,1,1,yind2[1]), c(-1,-1,1,length(yind2)))
gph <- abind::abind(gph1,gph2)[,length(glat):1,];rm(gph1,gph2)
gtimesel <- gtime[c(yind1, yind2)]
nc_close(ncg)

ncg_anom <- nc_open(gph_anom_file)
gtime <- as.Date(ncvar_get(ncg_anom,"time")/24, origin = "1800-01-01")
glon <- ncvar_get(ncg_anom,"longitude")
glon[glon > glon[length(glon)]] <- glon[glon > glon[length(glon)]] -360
glat <- ncvar_get(ncg_anom,"latitude");glat <- rev(glat)
yind1 <- which(year(gtime) == 1862);yind2 <- which(year(gtime) == 2011)
gph1 <- ncvar_get(ncg_anom,"gh", c(1,1,1,yind1[1]), c(-1,-1,1,length(yind1)))
gph2 <- ncvar_get(ncg_anom,"gh", c(1,1,1,yind2[1]), c(-1,-1,1,length(yind2)))
gph_anom <- abind::abind(gph1,gph2)[,length(glat):1,];rm(gph1,gph2)
nc_close(ncg_anom)

nct <- nc_open(t2m_file)
ttime <- as.Date(ncvar_get(nct,"time")/24, origin = "1800-01-01")
tlon <- ncvar_get(nct,"longitude")
tlon[tlon > tlon[length(tlon)]] <- tlon[tlon > tlon[length(tlon)]] -360
tlat <- ncvar_get(nct,"latitude");tlat <- rev(tlat)
yind1 <- which(year(ttime) == 1862);yind2 <- which(year(ttime) == 2011)
t2m1 <- ncvar_get(nct,"t2m", c(1,1,yind1[1]), c(-1,-1,length(yind1)))
t2m2 <- ncvar_get(nct,"t2m", c(1,1,yind2[1]), c(-1,-1,length(yind2)))
t2m <- abind::abind(t2m1,t2m2)[,length(tlat):1,];rm(t2m1,t2m2)
timesel <- ttime[c(yind1, yind2)]
nc_close(nct)

##### plot atmos figure  #######
brks <- seq(-200,200,25)
cols <- rev(colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(length(brks)-1))
marscale <- c(1,0.5,1.5,2); cx = 1


##### figure for the two springs
png(paste0(plotdir,"Figure_Warmsprings_Atmos_",Sys.Date(), ".png"), width = 2900, height = 1500, res = 300, pointsize = 8)
layout(matrix(1:8, nrow = 2, ncol = 4, byrow = T), heights = c(1,1), widths = c(rep(1,3),0.18))
par(mar = c(0.5,0.5,1,0.5), oma = c(2,2,1,1.5))

tind <- year(timesel) == 1862 & month(timesel) %in% c(3:5)
image(glon,glat,apply(gph_anom[,,tind],1:2,mean), xlim = c(-20,50), breaks = brks, col = cols, xaxt= "n", yaxt = "n")
maps::map("world",add = T, col = "gray75")
maps::map("world",region = "Switzerland", add = T, col = "red")
contour(glon,glat, apply(gph[,,tind],1:2,mean), add = T, col = "gray25")
mtext(side = 3, text = "MAM 1862")

spells <- list(c("1862-04-13","1862-04-14","1862-04-15","1862-04-16"),
               c("1862-04-26","1862-04-27","1862-04-28","1862-04-29","1862-04-30","1862-05-01","1862-05-02","1862-05-03","1862-05-04"))

for(day in 1:length(spells)){
  mat <- apply(gph_anom[,,as.character(timesel) %in% spells[[day]]],1:2,mean)
  mat[mat < brks[1]] <- brks[1];mat[mat > brks[length(brks)]] <- brks[length(brks)]
  image(glon,glat,mat, xlim = c(-20,50), breaks = brks, col = cols, xaxt= "n", yaxt = "n")
  maps::map("world",add = T, col = "gray75")
  maps::map("world",region = "Switzerland", add = T, col = "red")
  contour(glon,glat, apply(gph[,,as.character(timesel) %in% spells[[day]]],1:2,mean), add = T, col = "gray25")
  mtext(side = 3, text = paste0(spells[[day]][1], " to ",substr(spells[[day]][length(spells[[day]])],6,10)))
}

image_scale(breaks=brks,col=cols, axis.pos=4, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=4,ats=brks[seq(1,length(brks),2)],labs=brks[seq(1,length(brks),2)],
                            cex=cx,axis.pos2=4,line=3, adj=-0.07))

par(mar = c(0.5,0.5,1,0.5))
for(mm in 3:5){
  tind <- year(timesel) == 2011 & month(timesel) %in% mm
  image(glon,glat,apply(gph_anom[,,tind],1:2,mean), xlim = c(-20,50), breaks = brks, col = cols, xaxt= "n", yaxt = "n")
  maps::map("world",add = T, col = "gray75")
  maps::map("world",region = "Switzerland", add = T, col = "red")
  contour(glon,glat, apply(gph[,,tind],1:2,mean), add = T, col = "gray25")
  mtext(side = 3, text = paste0(c("March","April","May")[mm-2]," 2011"))
}

image_scale(breaks=brks,col=cols, axis.pos=4, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
            useraxis = list(axis.pos=4,ats=brks[seq(1,length(brks),2)],labs=brks[seq(1,length(brks),2)],
                            cex=cx,axis.pos2=4,line=3, adj=-0.07))

dev.off()
