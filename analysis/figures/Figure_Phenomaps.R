########  PLOT OF ALL PHENOLOGICAL INDICES #####
rm(list=ls())

### packages
library(ncdf4);library(lubridate)
library(viridis);library(RColorBrewer)

### functions
source("R/helpfuns.R")
source("R/image_funs.R")

### directories
plotdir <- "manuscript/figures"

indices <- list(maesh13d = list(file ="data/04_pheno_recon/CH_M1_maesh13d_1763-01-01-2020-12-31_2023-08-02.nc", var = "DOY", orig = "1763-01-01"),
                mprua65d = list(file ="data/04_pheno_recon/CH_PTTs_mprua65d_1763-01-01-2020-12-31_2023-08-02.nc", var = "DOY", orig = "1763-01-01"),
                mfags13d = list(file ="data/04_pheno_recon/CH_M1_mfags13d_1763-01-01-2020-12-31_2023-08-02.nc", var = "DOY", orig = "1763-01-01"))

nams <- names(indices)

### read topo
nc <- nc_open("../swiss_arm/code/data_input/topo.swiss02_ch01r.swisscors.nc")
topo <- ncvar_get(nc, "height")
topo_help <- topo
topo_help[topo < 1800] <- NA
topo_help[topo >= 1800] <- 1
topo_N <-ncvar_get(nc, "chy") + 1000000
topo_E <-ncvar_get(nc, "chx") + 2000000
nc_close(nc)

for(ii in 1:length(indices)){
  (var <- names(indices)[ii])
  nc <-nc_open(indices[[ii]]$file)
  time <- as.Date(ncvar_get(nc,"time"), origin = indices[[ii]]$orig)
  dat <- ncvar_get(nc,indices[[ii]]$var)
  assign(paste0(var,"_list"),list(dat = dat, time = time))
}


## calc climatology for each period
periods <- list(per1 = c(1781,1810), per2 = c(1811,1840), per3 = c(1841,1870) ,per_ind = c(1871,1900),
                per5 = c(1901,1930), per6 = c(1931,1960), per7 = c(1961,1990), per8 = c(1991,2020))

na_th <- function(x){sum(is.na(x))/length(x) > 0.8}
calc <- F

if(calc){
  for (ind in 1:length(indices)){
    phenolist <- get(paste0(nams[ind],"_list"))
    climmat <- array(sapply(periods, function(x){
      xind <- which(year(phenolist$time) %in% c(x[1]:x[2]))
      return(apply(phenolist$dat[,,xind],1:2,mean, na.rm = T))
    }), dim = c(370,240,length(periods)))
    assign(paste0("clim_",nams[ind]), climmat)
    save(climmat, file = paste0("data/04_pheno_recon/clim_",nams[ind],".RData"))

    ## calculate t-test
    pmat <- array(sapply(periods, function(p){
      print(p)
      xind <- which(year(phenolist$time) %in% c(p[1]:p[2]))
      xind_pre <- which(year(phenolist$time) %in% c(periods$per_ind[1]:periods$per_ind[2]))
      xdat <- plyr::alply(phenolist$dat[,,xind],c(1,2))
      ydat <- plyr::alply(phenolist$dat[,,xind_pre],c(1,2))
      mat <- mapply(function(x,y){if (any(na_th(x),na_th(y))){return(NA)} else{t.test(x,y, na.action ="na.pass")$p.value} },x = xdat, y = ydat)
      return(array(data = mat, dim = c(370,240)))
    }), dim = c(370,240,length(periods)))
    assign(paste0("pval_",nams[ind]), pmat)
    save(pmat, file = paste0("data/04_pheno_recon/pval_",nams[ind],".RData"))
  }

} else{

  for(ind in 1:length(indices)){
    load(paste0("data/04_pheno_recon/clim_",nams[ind],".RData"));assign(paste0("clim_",nams[ind]), climmat)
    load(paste0("data/04_pheno_recon/pval_",nams[ind],".RData"));assign(paste0("pval_",nams[ind]), pmat)
  }
}

rm(climmat)

### temperature indices
brks <- seq(85,150,5)
cols <- colorRampPalette(c("#0E6134","#2BAC66","#ABCF63","#E8F69E","#D7E3EE","#8FB3FF","#7F97FF","#0570B0","#023858"))(length(brks)-1)

marscale <- c(2,1,2,5);cx <- 1.5;adj = 0.5;cx1 = 1.5

days <- substr(seq(as.Date("1900-01-01"),as.Date("1900-07-01"), by = "day"),6,10)
labs <- days[brks[seq(1,length(brks),2)]]

png(paste0(plotdir,"/Figure_Phenomaps_",Sys.Date(),".png"), width =  length(periods) * 700, height = length(indices) * 550, res = 300, pointsize = 7)
layout(matrix(1:((length(periods) + 1) *3), nrow = 3, ncol = length(periods)+1, byrow = T), widths = c(rep(1,length(periods)),0.26), heights = 1)
par(oma = c(1,2,2,1))
for(pv in 1:3){
  par(mar = c(1,1,1,1))
  dat<- get(paste0("clim_",nams[pv]))
  dat[dat > brks[length(brks)]] <- brks[length(brks)]

  for(i in 1:length(periods)) {
    image(1:370,1:240,dat[,,i], xaxt  = "n", yaxt = "n", bty = "n", breaks = brks, col = cols)
    image(1:370,1:240,topo_help, col = "gray65", add = T)
    if(i == 1) mtext(side = 2, nams[pv], cex = cx1)
    if(i == 1) mtext(side = 3, paste0(letters[pv],")"), cex = cx1, adj = 0.01)
    if(pv == 1) mtext(side = 3, text = paste0(periods[[i]][1]," - ", periods[[i]][2]), cex = cx1, adj = 0.5)
  }
  image_scale(breaks=brks,col=cols, axis.pos=4, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
              useraxis = list(axis.pos=4,ats=brks[seq(1,length(brks),2)],labs=labs, cex=cx,axis.pos2=1,line=3, adj=-0.07))

}
dev.off()


### Anomalies plot
brks <- seq(-30,30,5)
cols <- colorRampPalette(rev(c("#6F5741","#B48222","white","#619E40","#326329")))(length(brks)-1)

png(paste0(plotdir,"/Figure_Phenomaps_Anomalies_",Sys.Date(),".png"), width =  length(periods) * 700, height = length(indices) * 550, res = 300, pointsize = 7)
layout(matrix(1:((length(periods) + 1) *3), nrow = 3, ncol = length(periods)+1, byrow = T), widths = c(rep(1,length(periods)),0.26), heights = 1)
par(oma = c(1,2,2,1))

for(pv in 1:3){
  par(mar = c(1,1,1,1))
  dat<- sweep(get(paste0("clim_",nams[pv])), MARGIN = c(1,2), STATS = get(paste0("clim_",nams[pv]))[,,which(names(periods)=="per_ind")])
  dat[dat > brks[length(brks)]] <- brks[length(brks)]
  dat[dat < brks[1]] <- brks[1]

  for(i in 1:length(periods)) {
    image(1:370,1:240,dat[,,i], xaxt  = "n", yaxt = "n", bty = "n", breaks = brks, col = cols)
    if(i == 1) mtext(side = 2, nams[pv], cex = cx1)
    if(i == 1) mtext(side = 3, paste0(letters[pv],")"), cex = cx1, adj = 0.01)
    if(pv == 1) mtext(side = 3, text = paste0(periods[[i]][1]," - ", periods[[i]][2]), cex = cx1, adj = 0.5)
  }
  image_scale(breaks=brks,col=cols, axis.pos=4, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
              useraxis = list(axis.pos=4,ats=brks[seq(1,length(brks),2)],labs=brks[seq(1,length(brks),2)], cex=cx,axis.pos2=1,line=3, adj=-0.07))

}
dev.off()



####### plot pvalues of differnces #######
brks_sig <- c(0,0.05,seq(0.1,1, 0.1))
cols_sig <- colorRampPalette(c("#225ea8", "#41b6c4", "#a1dab4", "#ffffcc", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026"))(length(brks_sig)-1)

png(paste0(plotdir,"/Figure_Phenomaps_Pvalues_",Sys.Date(),".png"), width =  length(periods) * 700, height = length(indices) * 550, res = 300, pointsize = 7)
layout(matrix(1:((length(periods) + 1) *3), nrow = 3, ncol = length(periods)+1, byrow = T), widths = c(rep(1,length(periods)),0.26), heights = 1)
#layout(matrix(1:((length(periods) + 1) *2), nrow = 2, ncol = length(periods)+1, byrow = T), widths = c(rep(1,length(periods)),0.26), heights = 1)
par(oma = c(1,2,2,1))
for(pv in 1:3){
  par(mar = c(1,1,1,1))
  dat<- get(paste0("pval_",nams[pv]))
  dat[dat > brks_sig[length(brks_sig)]] <- brks_sig[length(brks_sig)]

  for(i in 1:length(periods)) {
    image(1:370,1:240,dat[,,i], xaxt  = "n", yaxt = "n", bty = "n", breaks = brks_sig, col = cols_sig)
    image(1:370,1:240,topo_help, col = "gray65", add = T)
    if(i == 1) mtext(side = 2, nams[pv], cex = cx1)
    if(i == 1) mtext(side = 3, paste0(letters[pv],")"), cex = cx1, adj = 0.01)
    if(pv == 1) mtext(side = 3, text = paste0(periods[[i]][1]," - ", periods[[i]][2]), cex = cx1, adj = 0.5)
  }
  image_scale(breaks=brks_sig,col=cols_sig, axis.pos=4, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
              useraxis = list(axis.pos=4,ats=brks_sig[seq(1,length(brks_sig),2)],labs=brks_sig[seq(1,length(brks_sig),2)], cex=cx,axis.pos2=1,line=3, adj=-0.07))
}
dev.off()







