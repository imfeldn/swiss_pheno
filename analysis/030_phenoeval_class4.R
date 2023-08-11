##### COMPARE MODELS TO INDEPENDENT CLASS 4 OBSERVATIONS ########
### this are stations that have not been used in the calibrations

rm(list=ls())

### packages
library(ncdf4);library(mmand);library(MASS)

### functions
source("R/helpfuns.R")

### directories
indir <- "data/03_calibrations/"
plotdir <- "manuscript/supplementary/"

### input
subDir <- "sa_run_2023-08-10_class123_lateperiod"
phenovals <- c("mprua65d","mfags13d","maesh13d")#,"mtaro65d","manen65d"
allyrs <- 1951:2020

#pv <- "mfags13d"

for(pv in phenovals){

  ## load pheno & meteo and create phenor object
  (load(paste0("data/02_pheno_net/pheno_meta_",pv,"_class4.RData")))
  (load(paste0("data/02_pheno_net/meteo_",pv,"_class4.RData")))

  stns_pv <- unique(pheno_meta$nat_abbr)
  stn_list_class4 <- list()

  # convert data to phenor format
  for(stn in 1:length(stns_pv)){

    id <- stns_pv[stn]
    ind <- which(pheno_meta$nat_abbr == id & pheno_meta$reference_year < 2021)
    locs <- as.numeric(unlist(pheno_meta[ind[1],c("lat","lon")]))
    stn_list_class4[[stn]]  <- list(site = id, location = locs,  doy = c(-102:-1,1:263),
                             ltm = rep(15,365), # I'm not sure what this is for
                             transition_dates = pheno_meta$doy[ind],
                             year = pheno_meta$reference_year[ind],
                             Ti = t(meteo_dat[rownames(meteo_dat) == id,allyrs %in% pheno_meta$reference_year[ind],]),
                             Tmini = t(meteo_dat[rownames(meteo_dat) == id,allyrs %in% pheno_meta$reference_year[ind],]) - 3, # this is artificial data for tmin, but never used in the models
                             Tmaxi = t(meteo_dat[rownames(meteo_dat) == id,allyrs %in% pheno_meta$reference_year[ind],]) + 3, # this is artificial data for tmax, but never used in the models
                             Li = matrix(rep(daylength(c(264:365,1:263),locs[1]), length(ind)), nrow = 365, ncol = length(ind)))


  }

  names(stn_list_class4) <- stns_pv


  ## load models
  (load(paste0("data/03_calibrations/",subDir ,"/",pv,"/model_comparison_",pv,".RData")))
  mods <- names(mod_out$modelled)
  seeds <- 1:nrow(mod_out$modelled$LIN$parameters)

  ## predict class 4 stations
  eval_mods_class4 <- array(NA, dim=c(length(stn_list_class4),length(mods),4,length(seeds)))
  dimnames(eval_mods_class4) <- list(names(stn_list_class4),mods, c("rmse","bias","pcorr","nse"),seeds)

  for (mod in 1:length(mod_out$modelled)){

    for(ss in 1:length(seeds)){

      pred <- pr_predict(par = mod_out$modelled[[mod]]$parameters[ss,], data = stn_list_class4, model = names(mod_out$modelled)[mod])
      pred[pred > 1000] <- NA
      t <- 0

      for(stn in 1:length(stn_list_class4)){
        stnind <- (t + 1):(t + length(stn_list_class4[[stn]]$transition_dates))
        t <- stnind[length(stnind)]
        obs <- stn_list_class4[[stn]]$transition_dates
        eval_mods_class4[stn,mod,,ss] <- round(eval_recon(rec = pred[stnind], obs)[c("rmse","bias","pcorr","nse"),],3)
      }
    }
  }
  save(eval_mods_class4, file = paste0("data/03_calibrations/",subDir,"/",pv,"/eval_mods_class4_",pv,".RData"))



  # plot the evaluation ###
  print("plot eval")

  nrclass4 <- nrow(eval_mods_class4)
  nc <- length(seeds) * length(mods)
  cols <- rep(brewer.pal(n = length(mods), "YlGnBu"), each = length(seeds))
  class4 <- aperm(eval_mods_class4, perm = c(1,4,3,2))

  lwd <- 0.5
  png(paste0("data/03_calibrations/",subDir,"/",pv,"/class4_",pv,"_",Sys.Date(),".png"), width = 3200, height = 900, res = 300, pointsize = 6)
  par(mfcol = c(1,4), mar = c(1,2,0,1), oma = c(2,2,2,2))
  rr <- range(class4[,,1,],na.rm = T)
  if (rr[2] > 25) rr[2] <- 25
  df <- data.frame(matrix(class4[,,1,],nrow = nrclass4, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  mtext(side = 3,line = 0.1, "RMSE")
  abline( h = seq(0,25,5), lty = 3, col = "gray75", lwd = lwd)
  axis(side = 1, at = seq(2,nc,length(seeds)), mods)

  rr <- range(class4[,,2,],na.rm = T)
  df <- data.frame(matrix(class4[,,2,],nrow = nrclass4, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  mtext(side = 3,line = 0.1, "Mean bias")
  abline( h = seq(-25,25,5), lty = 3, col = "gray75", lwd = lwd)
   axis(side = 1, at = seq(2,nc,length(seeds)), mods)

   rr <- range(class4[,,3,],na.rm = T)
   df <- data.frame(matrix(class4[,,3,],nrow = nrclass4, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  mtext(side = 3,line = 0.1, "Pearson correlation")
  abline( h = seq(-1,1,0.2), lty = 3, col = "gray75", lwd = lwd)
  axis(side = 1, at = seq(2,nc,length(seeds)), mods)

  rr <- range(class4[,,4,],na.rm = T)
  if (rr[1] < -3) rr[1] <- -3
  df <- data.frame(matrix(class4[,,4,],nrow = nrclass4, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  mtext(side = 3,line = 0.1, "NSE")
  abline( h = seq(-4,1,0.5), lty = 3, col = "gray75", lwd = lwd)
  axis(side = 1, at = seq(2,nc,length(seeds)), mods)
  dev.off()

  # for(ss in 1:length(stations)){
  #
  #   sscoord <- pheno_meta[which(pheno_meta$nat_abbr == stations[ss])[1],c("E","N")]
  #
  #   stnind <- c(which.min(abs(E - sscoord[1,1])),which.min(abs(N - sscoord[1,2])))
  #   stn_recon <- hcn[stnind[1],stnind[2],]
  #
  #   stn_obs <- pheno_meta[pheno_meta$nat_abbr == stations[ss],c("reference_year","doy")]
  #
  #   met_lti <- c()
  #   met_lti[1] <- round(cor(stn_recon[time %in% stn_obs$reference_year], stn_obs$doy[stn_obs$reference_year %in% time], use = "na.or.complete"),2)
  #   met_lti[2] <- round(mean((stn_recon[time %in% stn_obs$reference_year] - stn_obs$doy[stn_obs$reference_year %in% time])^2, na.rm = T)^0.5,2)
  #
  #   png(paste0(plotdir,"eval_",names(ids)[pv],"_",stations[ss],"_",Sys.Date(),".png"), width = 600, height = 400)
  #   par(oma = c(3,1,3,1), cex = 1.2, mar = c(1.5,3,1,1))
  #   days <- seq(as.Date("1900-01-01"),as.Date("1900-07-01"), by = "day")
  #   rr <- range(c(stn_obs$doy,stn_recon), na.rm = T)
  #   plot(stn_obs$reference_year, stn_obs$doy, type = "l", ylab = "", xlab = "", yaxt = "n", xlim = c(1950, 2020),
  #        ylim = rr + c(-5,5), col = t_col("plum",240), lwd = 2)
  #   points(stn_obs$reference_year, stn_obs$doy, type = "p", ylab = "", xlab = "", yaxt = "n", col = t_col("plum",180), pch = 19)
  #   points(1763:2020, stn_recon, type = "l", col =t_col("darkblue",120), lwd = 2)
  #   legend("topright", legend = c("obs","rec"), col = c("plum", "darkblue"), lwd = 2, bty = "n")
  #   abline(h = which(substr(days,9,10) == "01"), lty = 2, col = "gray75")
  #   text(1962, y = rr[2]-2, paste0("PCor: ", met_lti[c(1)]))
  #   text(1962, y = rr[2]-6, paste0("RMSE: ",  met_lti[c(2)]))
  #   ind <- substr(days,9,10) %in% c("01","15")
  #   axis(side = 2, at = which(ind) ,ind , las = 2)
  #   mtext(side = 3, paste0(names(ids)[pv], stations[ss]), line  = -1.2, cex = 1.2)
  #   dev.off()
  #
  # }
}

