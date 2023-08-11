##### PLOT LATE SPRING Conditions #######
rm(list=ls())

### packages
library(ncdf4);library(lubridate);library(RColorBrewer)
library(XLConnect)

### helpfer functions
source("R/image_funs.R")
source("R/plot_winds.R")
#source("../../ncode/R/randomfuns.R")

### directores
plotdir <- "manuscript/figures/"

### adjust spring index still ###
(load("../swiss_indices/5_output/4_pheno_recon/Springindex_phenorecon_2023-06-13.RData"))

early <- df_springs_yr$years[order(df_springs_yr$ranks)][1:10]
late <- rev(df_springs_yr$years[order(df_springs_yr$ranks)])[2:11]

latesprgs <- late#c(1785,1837,1853)

### load stormniess
(load("analysis/storminess_MAM_lanzcos_ref19611990_2023-08-10.RData"))
### read inventory
inventory <- readWorksheetFromFile("../swiss_arm/code/data_input/arm_stations_inventory_withnbcn.xlsx",
                                   sheet=1, colTypes="character",startRow=2)

inv <- inventory[inventory$ID %in% colnames(sd_anom_MAM),c("ID","lon_hist","lat_hist")]
inv <- inv[match(colnames(sd_anom_MAM), inv$ID),]

######## weather types during early spring months ########
load("../swiss_arm/code/data_input/analogue_stationdata.RData")
wts <- TOT[,1:15];rm(TOT)
wtnames <- c("NE","WSW","W","E","HP","N","WC")

wts_latest  <- wts[year(wts$date) %in% latesprgs & month(wts$date) %in% c(3:5),]

agg_wts_latest <- sapply(1:nrow(wts_latest), function(x){
  wttype <- which(cumsum(as.numeric(wts_latest[x,seq(3,15,2)])) > 0.9)[1]
  print(wttype)
  tab <- data.frame(matrix(NA,ncol = 8))
  tab[1:(wttype +1)] <- wts_latest[x, c(1,seq(2,15,2)[1:wttype])]
  return(tab)
})

agg_wts_latest <- data.frame(matrix(unlist(agg_wts_latest), nrow = nrow(wts_latest), ncol = 8, byrow = T))
colnames(agg_wts_latest) <- c("date",paste0("wt",1:7))
agg_wts_latest$date <- as.Date(agg_wts_latest$date, "1970-01-01")

for(ss in 1:length(latesprgs)){
  assign(paste0("wts_",latesprgs[ss]),table(unlist(agg_wts_latest[year(agg_wts_latest$date) == latesprgs[ss],2:8])))
  wts_perc <- round(get(paste0("wts_",latesprgs[ss]))/sum(get(paste0("wts_",latesprgs[ss]))) * 100,2)
  assign(paste0("wts_",latesprgs[ss],"perc"), wts_perc)
}


### all wts
wts_all <- wts[month(wts$date) %in% c(3:5),]

wts_all[is.na(wts_all)] <- 1
all_wts <- sapply(1:nrow(wts_all), function(x){
  wttype <- which(cumsum(as.numeric(wts_all[x,seq(3,15,2)])) > 0.9)[1]
  print(wttype)
  tab <- data.frame(matrix(NA,ncol = 8))
  tab[1:(wttype +1)] <- wts_all[x, c(1,seq(2,15,2)[1:wttype])]
  return(tab)
})

wts_all_all <- table(unlist(all_wts[2:8,]))
wts_all_all_perc <- wts_all_all/sum(wts_all_all) * 100

### late wts
wts_all_lates <- wts[month(wts$date) %in% c(3:5) & year(wts$date) %in% late,]

late_wts <- sapply(1:nrow(wts_all_lates), function(x){
  wttype <- which(cumsum(as.numeric(wts_all_lates[x,seq(3,15,2)])) > 0.9)[1]
  print(wttype)
  tab <- data.frame(matrix(NA,ncol = 8))
  tab[1:(wttype +1)] <- wts_all_lates[x, c(1,seq(2,15,2)[1:wttype])]
  return(tab)
})

wts_all_late <- table(unlist(late_wts[2:8,]))
wts_all_late_perc <- wts_all_late/sum(wts_all_late) * 100

###### plots #################

ylim <- c(43.8,50.4)
xlim <- c(5.4,12.5)

png(paste0(plotdir,"Figure_Latesprings_WTS_",Sys.Date(), ".png"), width = 1800,height = 1000, res = 300, pointsize = 6)
layout(matrix(c(1:8,9,9,9,9), nrow = 3, ncol = 4, byrow = T), width = c(1), heights = c(1,1,0.18))
par(mar = c(1,0,3,0), oma = c(2,3,1,2))
## weather types
barplot(wts_1785perc - wts_all_all_perc, names = wtnames, space = 0, ylim = c(-12,20));box()
mtext(side = 3, text = 1785, line = -1.5);abline(h = 0)
mtext(side = 3, text = "a)", adj = 0.01)
axis(side = 1, at = seq(1:7), labels = FALSE, tick = TRUE, line = NA, outer = FALSE,lwd.ticks = 1, col.ticks = "black")
barplot(wts_1837perc - wts_all_all_perc, names = wtnames, space = 0, ylim = c(-12,20),  yaxt = "n");box();mtext(side = 3, text = 1837, line = -1.5);abline(h = 0)
axis(side = 1, at = seq(1:7), labels = FALSE, tick = TRUE, line = NA, outer = FALSE,lwd.ticks = 1, col.ticks = "black")
barplot(wts_1853perc - wts_all_all_perc, names = wtnames, space = 0, ylim = c(-12,20), yaxt = "n");box();mtext(side = 3, text = 1853, line = -1.5);abline(h = 0)
axis(side = 1, at = seq(1:7), labels = FALSE, tick = TRUE, line = NA, outer = FALSE,lwd.ticks = 1, col.ticks = "black")
barplot(wts_all_late_perc  - wts_all_all_perc, names = wtnames, space = 0, ylim = c(-12,20), yaxt = "n");box();mtext(side = 3, text = "Late springs", line = -1.5);abline(h = 0)
axis(side = 1, at = seq(1:7), labels = FALSE, tick = TRUE, line = NA, outer = FALSE,lwd.ticks = 1, col.ticks = "black")
## storminess
brks_storm <- seq(-1,1,0.2)
cols_storm <- colorRampPalette(c("#800000", "#8C5000", "#C58917", "#FFA500", "#008000", "#307D7E", "#2B65EC", "#0000A0"))(length(brks_storm)-1)

for(ii in c(1785,1837,1853)){
  selyr <- which(1763:2020 == ii)
  image(seq(5,15,1),seq(41,55,1),matrix(NA, length(seq(5,15,1)), length(seq(41,55,1))),xlab="", ylab="",
        asp=1,xaxt="n",yaxt="n", xlim = xlim, ylim = ylim)
  maps::map("world", add=T, col="gray25")
  colp <- cols_storm[as.numeric(cut(round(sd_anom_MAM[selyr,],2),breaks = brks_storm))]
  nna <- which(!is.na(colp))
  points(inv$lon_hist[nna], inv$lat_hist[nna], pch = 15, cex = 2, col = colp[nna])
  mtext(side =3, line = -1.5, text = ii)
  if(ii == 1785) mtext(side = 3, text = "b)", adj = 0.01)
}

image(seq(5,15,1),seq(41,55,1),matrix(NA, length(seq(5,15,1)), length(seq(41,55,1))),xlab="", ylab="",
      asp=1,xaxt="n",yaxt="n", xlim = xlim, ylim = ylim)
maps::map("world", add=T, col="gray25")
colp <- cols_storm[as.numeric(cut(round(colMeans(sd_anom_MAM[1763:2020 %in% late[1:10],], na.rm = T),2),breaks = brks_storm))]
nna <- which(!is.na(colp))
points(inv$lon_hist[nna], inv$lat_hist[nna], pch = 15, cex = 2, col = colp[nna])
mtext(side =3, line = -1.5, text = "Late springs")

par(mar=c(0,0,0,0))
image(seq(5,15,1),seq(41,55,1),matrix(NA, length(seq(5,15,1)), length(seq(41,55,1))),xlab="", ylab="",
      asp=1,xaxt="n",yaxt="n", bty = "n")
legend("top",legend = paste0(brks_storm[1:(length(brks_storm)-1)]," to ",brks_storm[2:length(brks_storm)]) ,
       col = cols_storm, pch = 15, bty = "n", ncol = 5, pt.cex = 2, cex = 1.5)
dev.off()
