##### PLOT LATE SPRING INDICES #######
rm(list=ls())

### packages
library(ncdf4);library(lubridate);library(RColorBrewer)

### helpfer functions
source("R/image_funs.R")

### directores
plotdir <- "manuscript/figures/"

periods <- list(per1 = c(1781,1810),per2 = c(1841,1870),per3 = c(1991,2020))

altth <- 1600

nc <- nc_open("../swiss_arm/code/data_input/topo.swiss02_ch01r.swisscors.nc")
topo <- ncvar_get(nc, "height")
topo_help <- topo
topo_help[topo >= altth] <- 10
topo_help[topo < altth] <- 1
topo_N <-ncvar_get(nc, "chy") + 1000000
topo_E <-ncvar_get(nc, "chx") + 2000000
nc_close(nc)

#### get warm/cold years for the swiss plateau....
sfc <- c("Alp", "AlpsS","Jura", "Swiss Plateau", "Pre-Alps")
regind <- which(sfc == "Swiss Plateau")
nc <- nc_open("data/01_climate_indices/tempmean/fldmean_regions_TEMPMEAN_season_1763-2020.nc")
tmean_seas <- ncvar_get(nc, "temp")
tmean_seas_time <- as.Date(ncvar_get(nc, "time"), origin = "1970-01-01")

### create a table for comparison ####
df_temp_all <- data.frame(yr = 1763:2020, temp = tmean_seas[regind,month(tmean_seas_time) == 4])
df_temp_all_present <- df_temp_all[df_temp_all$yr %in% periods[[3]][1]:periods[[3]][2],]
df_temp_all_preindust <- df_temp_all[df_temp_all$yr %in% periods[[2]][1]:periods[[2]][2],]
df_temp_all_hist <- df_temp_all[df_temp_all$yr %in% periods[[1]][1]:periods[[1]][2],]

ord_all <- order(df_temp_all$temp)
View(df_temp_all[ord_all,])

ord_pres <- order(df_temp_all_present$temp)
df_temp_all_present[ord_pres,]
View(data.frame(yr = df_temp_all_present[ord_pres,1], anom = round(df_temp_all_present[ord_pres,2] - mean(df_temp_all_present[ord_pres,2] ),1)))

ord_preind <- order(df_temp_all_preindust$temp)
df_temp_all_preindust[ord_preind,]
View(data.frame(yr = df_temp_all_preindust[ord_preind,1], anom = round(df_temp_all_preindust[ord_preind,2] - mean(df_temp_all_preindust[ord_preind,2] ),2)))

ord_hist <- order(df_temp_all_hist$temp)
df_temp_all_hist[ord_hist,]
View(data.frame(yr = df_temp_all_hist[ord_hist,1], anom = round(df_temp_all_hist[ord_hist,2] - mean(df_temp_all_hist[ord_hist,2] ),1)))

preind <- rev(df_temp_all_preindust[ord_preind,"yr"])[1:3]
present <- rev(df_temp_all_present[ord_pres,"yr"])[1:3]
hist <- rev(df_temp_all_hist[ord_hist,"yr"])[1:3]


(load("data/01_climate_indices/swissplateau_stations_indices_2023-08-10.RData"))
mam_sw <- data.frame(index_list$tempmean[,c("year","MAM")])
View(mam_sw[order(mam_sw[,2]),])

mam_sw$anom4170 <- mam_sw[,2] - mean(mam_sw[mam_sw$year %in% periods$per2[1]:periods$per2[2],2])
mam_sw$anom9120 <- mam_sw[,2] - mean(mam_sw[mam_sw$year %in% periods$per3[1]:periods$per3[2],2])

### load indices
springs <- c(1862,2011)

indices <- list(tempmean = list(file = "data/01_climate_indices/tempmean/TEMPMEAN_season_1763-2020.nc", var = "temp", orig = "1970-01-01", longnams = "Mean temperature"),
                mprua65d = list(file ="5_output/04_pheno_recon/CH_PTTs_mprua65d_1763-01-01-2020-12-31_2023-08-02.nc", var = "DOY", orig = "1763-01-01", longnams = "Cherry flowering"),
                wsdi = list(file ="data/01_climate_indices/hwdi/WSDI_season_1763-2020_2023-07-26.nc", var = "days", orig = "1970-01-01", longnams = "Warm spell index"),
                wetdays = list(file ="data/01_climate_indices/wd/WD_season_1763-2020.nc", var = "precip", orig = "1970-01-01", longnams = "Wet days"),
                snowdays = list(file ="data/01_climate_indices/snowdays/SNOWDAYS_season_1763-2020.nc", var = "temp", orig = "1970-01-01", longnams = "Snowfall days"))

for(ii in 1:length(indices)){
  (var <- names(indices)[ii])
  nc <- nc_open(indices[[ii]]$file)
  time <- as.Date(ncvar_get(nc,"time"), origin = indices[[ii]]$orig)
  tind <- which(year(time) %in% springs )
  if(length(tind) > length(springs)){tind <- which(year(time) %in% springs  & month(time) %in% c(3:5))}
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
  tind <- switch((indices[[ii]]$var == "DOY") + 1, which(year(time) %in% 1871:1900  & month(time) %in% c(3:5)),which(year(time) %in% 1871:1900))
  assign(paste0(var,"_mat_clim"),apply(ncvar_get(nc,indices[[ii]]$var)[,,tind],1:2,mean))
  nc_close(nc)
}


brksp_list <- list(tempmean = seq(-5,5,0.5),mprua = seq(-24,24,4),hwdi = seq(-24,24,4), wetdays = seq(-24,24,4), snowdays = seq(-24,24,4))
cols_list <- list(rev(colorRampPalette(brewer.pal(n =5, name = "RdBu"))(length(brksp_list[[1]])-1)),
                  colorRampPalette(c("#0E6134","#2BAC66","#ABCF63","#E8F69E","#D7E3EE","#8FB3FF","#7F97FF","#0570B0","#023858"))(length(brksp_list[[2]])-1),
                  rev(colorRampPalette(brewer.pal(n =5, name = "RdBu"))(length(brksp_list[[3]])-1)),
                  colorRampPalette(c("peru","gray98","cornflowerblue"))(length(brksp_list[[4]])-1),
                  colorRampPalette(c("peru","gray98","cornflowerblue"))(length(brksp_list[[4]])-1))

lty <- 1; adj = 0.01
cx = 1; marscale = c(3,2,1,2); mline = 2.5


png(paste0(plotdir,"Figure_Warmsprings_",Sys.Date(),".png"), width = 3200, height = 1300, res = 300, pointsize = 7)
layout(matrix(1:(3*5), nrow = 3, ncol = 5, byrow = F), widths = 1, heights = c(1,1,0.22))
par(oma=c(1,1,1,2))

for(ii in c(1,3,4,5,2)){
  (ind <- names(indices)[ii])
  brksp <-  brksp_list[[ii]]
  colp <-  cols_list[[ii]]

  for(l in 1:length(springs)){
    par(mar=c(1,1,1,1))
    mat <- get(paste0(ind,"_mat"))[,,l] - get(paste0(ind,"_mat_clim"))
    mat[mat < brksp[1]] <- brksp[1]
    mat[mat > brksp[length(brksp)]] <- brksp[length(brksp)]

    if(ii == 2){
      mat[topo > altth] <- NA
      image(1:nrow(mat),1:ncol(mat),topo_help, bty = "n", yaxt = "n", xaxt = "n", col = c("gray90","gray45"), breaks = c(0,1,10), xlab = "", ylab = "")
      image(1:nrow(mat),1:ncol(mat),mat, bty = "n", yaxt = "n", xaxt = "n", col = colp, breaks = brksp, add = T)
    } else{
      image(1:nrow(mat),1:ncol(mat),mat, bty = "n", yaxt = "n", xaxt = "n", col = colp, breaks = brksp)
    }
    if(ii == 1) mtext(side = 2, text = springs[l])
    if(l == 1) mtext(side = 3, text = indices[[ii]]$longnams)
  }
  image_scale(breaks=brksp,col=colp, axis.pos=1, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
              useraxis = list(axis.pos=1,ats=brksp[seq(1,length(brksp),2)],labs=brksp[seq(1,length(brksp),2)], cex=cx,axis.pos2=1,line=2, adj=0))
}

dev.off()


