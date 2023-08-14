####### MAP PLOT EVALUATION STATION WISE #########


### packages
library(viridis);library(ncdf4)

### functions
source("R/image_funs.R")

### directories
indir <- "data/03_calibrations/"
plotdir <- "manuscript/supplementary/"

nc <- nc_open("../swiss_arm/code/data_input/topo.swiss02_ch01r.swisscors.nc")
topo <- ncvar_get(nc, "height")
topo_N <-ncvar_get(nc, "chy") + 1000000
topo_E <-ncvar_get(nc, "chx") + 2000000
nc_close(nc)

phenovals <- c("mprua65d","mfags13d","maesh13d")#,"mtaro65d","manen65d"
annealing_version <- "sa_run_2023-08-11_class123_longseries_lateperiod"

pv <- "maesh13d"

### read evaluation data
for(pv in phenovals){

  (load(paste0(indir,annealing_version,"/",pv,"/eval_mods_test_",pv,".RData")))
  assign(paste0(pv,"_test"), eval_mods_test)

  (load(paste0(indir,annealing_version,"/",pv,"/eval_mods_train_",pv,".RData")))
  assign(paste0(pv,"_train"), eval_mods_train)

  (load(paste0("data/02_pheno_net/stations_",pv,"_class123.RData")))
  assign(paste0(pv,"_meta"), station_list_sel)

  }

##### map plots of biases, etc. in space
brks_rms <- seq(0,25,2.5)
rmscols <-  viridis(length(brks_rms)-1)
brks_bias <- seq(-20,20,2.5)
biascols <- colorRampPalette(c("#00204DFF","white","darkred"))(length(brks_bias)-1)
brks_corr <- seq(-0.2,1,0.1)
corrcols <- rev(viridis(length(brks_corr)-1))
brks_nse <- seq(-3,1,0.2)
nsecols <-rev(viridis(length(brks_nse)-1))

lonvals <- c(7,9.3); latvals <- c(46.2,47.3)
marscale <- c(1,1,1,1)

met <- c("rms","bias","corr","nse")
mods <-  c("LIN", "TT", "TTs", "PTT", "PTTs","M1","M1s","AT")
cx <- 1

for(pv in phenovals){

  test <- get(paste0(pv,"_test"))
  train <- get(paste0(pv,"_train"))
  locs <- get(paste0(pv,"_meta"))

  png(paste0(plotdir, "evaluation_train_",annealing_version,"_",pv,"_",Sys.Date(), ".png"), height = 1600, width = 4600, res = 200, pointsize = 7)
  layout(matrix(1:(9*4), ncol=9, nrow=4, byrow=T),widths = c(rep(1,8),0.18), heights = c(1))
  par(oma = c(1,2,2,4), cex = 1.5)

  for(mm in 1:length(met)){

    brks <- get(paste0("brks_",met[mm]))
    cols <- get(paste0(met[mm],"cols"))

    dat_train <- train[,,mm,2]
    dat_train[dat_train <= brks[1]] <- brks[1] + 0.001
    dat_train[dat_train >= brks[length(brks)]] <- brks[length(brks)] - 0.001

    for(mod in 1:length(mods)){
      par(mar = c(0.2,0.2,0.2,0.2))
      image(topo_E, topo_N, topo, ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n", col = "gray")
      maps::map("world", add = T)
      colp = cols[as.numeric(cut(dat_train[,mod], breaks = brks))]
      points(locs[rownames(train),"KoordinatenE"], locs[rownames(train),"KoordinatenN"], pch = 21, bg = colp, cex = 1.5)
      if(mod == 1) mtext(side = 2, met[mm], cex = cx + 1, line = 1)
      if(mm == 1) mtext(side = 3, mods[mod], cex = cx + 1, line = 1)
       }

    image_scale(breaks=brks,col=cols, axis.pos=4, scale_lab="", equidist=T, add.axis=F,mar=marscale,
                useraxis = list(axis.pos=4,ats=1:length(brks),labs=brks, cex=1,axis.pos2=4,line=1.5, adj=1))

  }

  dev.off()


  png(paste0(plotdir, "evaluation_test_",annealing_version,"_",pv,"_",Sys.Date(), ".png"), height = 1600, width = 4600, res = 200, pointsize = 7)
  layout(matrix(1:(9*4), ncol=9, nrow=4, byrow=T),widths = c(rep(1,8),0.18), heights = c(1))
  par(oma = c(1,2,2,4), cex = 1.5)

  for(mm in 1:length(met)){

    brks <- get(paste0("brks_",met[mm]))
    cols <- get(paste0(met[mm],"cols"))

    dat_test <- test[,,mm,2]
    dat_test[dat_test <= brks[1]] <- brks[1] + 0.001
    dat_test[dat_test >= brks[length(brks)]] <- brks[length(brks)] - 0.001

    for(mod in 1:length(mods)){

      par(mar = c(0.2,0.2,0.2,0.2))
      image(topo_E, topo_N, topo, ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n", col = "gray")
      maps::map("world", add = T)
      colp = cols[as.numeric(cut(dat_test[,mod], breaks = brks))]
      points(locs[rownames(test),"KoordinatenE"], locs[rownames(test),"KoordinatenN"], pch = 21, bg = colp, cex = 1.5)
      if(mod == 1) mtext(side = 2, met[mm], cex = cx + 1, line = 0.1)
      if(mm == 1) mtext(side = 3, mods[mod], cex = cx + 1, line = 1)

    }

    image_scale(breaks=brks,col=cols, axis.pos=4, scale_lab="", equidist=T, add.axis=F,mar=marscale,
                useraxis = list(axis.pos=4,ats=1:length(brks),labs=brks, cex=1,axis.pos2=4,line=1.5, adj=1))

  }

  dev.off()


  png(paste0(plotdir, "evaluation_test_",annealing_version,"_",pv,"_",Sys.Date(), ".png"), height = 1600, width = 4600, res = 200, pointsize = 7)
  layout(matrix(1:(9*4), ncol=9, nrow=4, byrow=T),widths = c(rep(1,8),0.18), heights = c(1))
  par(oma = c(1,2,2,4), cex = 1.5)

  for(mm in 1:length(met)){

    brks <- get(paste0("brks_",met[mm]))
    cols <- get(paste0(met[mm],"cols"))

    dat_test <- test[,,mm,2]
    dat_test[dat_test <= brks[1]] <- brks[1] + 0.001
    dat_test[dat_test >= brks[length(brks)]] <- brks[length(brks)] - 0.001

    for(mod in 1:length(mods)){

      par(mar = c(0.2,0.2,0.2,0.2))
      image(topo_E, topo_N, topo, ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n", col = "gray")
      maps::map("world", add = T)
      colp = cols[as.numeric(cut(dat_test[,mod], breaks = brks))]
      points(locs[rownames(test),"KoordinatenE"], locs[rownames(test),"KoordinatenN"], pch = 21, bg = colp, cex = 1.5)
      if(mod == 1) mtext(side = 2, met[mm], cex = cx + 1, line = 0.1)
      if(mm == 1) mtext(side = 3, mods[mod], cex = cx + 1, line = 1)

    }

    image_scale(breaks=brks,col=cols, axis.pos=4, scale_lab="", equidist=T, add.axis=F,mar=marscale,
                useraxis = list(axis.pos=4,ats=1:length(brks),labs=brks, cex=1,axis.pos2=4,line=1.5, adj=1))

  }

  dev.off()

}


