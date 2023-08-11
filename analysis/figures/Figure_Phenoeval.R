##### PLOT PHENOLOGY EVALUATION #######
rm(list=ls())

### packages
library(ncdf4);library(lubridate);library(RColorBrewer)

### helpfer functions
source("R/image_funs.R")
source("R/plot_winds.R")
source("R/gaussian_wrap.R")

### directores
plotdir <- "manuscript/figures/"
annealing_version <- "sa_run_2023-08-09_class123"

### read topo in order to only use low elevation variables
nc <- nc_open("../swiss_arm/code/data_input/topo.swiss02_ch01r.swisscors.nc")
topo <- ncvar_get(nc, "height")
topobl <- which(topo < 1600)
topo_N <-ncvar_get(nc, "chy") + 1000000
topo_E <-ncvar_get(nc, "chx") + 2000000
nc_close(nc)

## pheno obs
(load("../swiss_indices/2_data/pheno_hist/eurohist_datlist.RData"))
liestal <- read.table("../swiss_indices/2_data/pheno_net/liestal.txt", header = T)
genf <- read.table("../swiss_indices/2_data/pheno_net/genf.txt", header = T)

path <- "5_output/4_pheno_recon/"
files_all <- list.files(path = path, pattern = c("CH_"))#maesh|mprua|mfags|mcora|mtaro|mcarp|manen|mlard|mtusf
files <- files_all[!grepl("swissplateau|index",files_all)]
pvs <- gsub("CH_M1*.|CH_PTTs*.|CH_PTT*.","",files)
pvs <- substr(pvs,1,8)

## coordinates from https://www.meteoschweiz.admin.ch/home/messwerte.html?param=messnetz-phaenologie&station=PGE
genf_coord <- c(2500470, 1116470)
liestal_coord <- c(2622208, 1260342)

ids = list(mprua65d = 662, maesh13d = 601)
pv = 1

## read recons
for(pv in 1:length(ids)){

  print(names(ids[pv]))
  ff <- files[grepl(names(ids[pv]),files)]
  nc <- nc_open(paste0(indir,"/",ff))
  hcn <- ncvar_get(nc,"DOY")
  print(str(hcn))
  hcn[,,1] <- NA
  time <- lubridate::year(as.Date(ncvar_get(nc,"time"), origin = "1763-01-01"))
  N <- ncvar_get(nc,"N")
  E <- ncvar_get(nc,"E")

  if(grepl("maesh",names(ids[pv]))){
    genfind <- c(which.min(abs(E - genf_coord[1])),which.min(abs(N - genf_coord[2])))
    stn_recon <- hcn[genfind[1],genfind[2],]
    stn_obs <- data.frame(reference_year = genf$Jahr, doy = genf$Tagnr)
    met_dat <- c()
    met_dat[1] <- round(cor(stn_recon[time %in% genf$Jahr], genf$Tagnr, use = "na.or.complete"),2)
    met_dat[2] <- round(mean((stn_recon[time %in% genf$Jahr] - genf$Tagnr)^2, na.rm = T)^0.5,2)
    met_dat[3] <- round(mean((stn_recon[time %in% genf$Jahr] - genf$Tagnr), na.rm = T),2)
    met_dat[4] <- round(cor(stn_recon[time %in% genf$Jahr[genf$Jahr < 1900]], genf$Tagnr[genf$Jahr < 1900], use = "na.or.complete"),2)
    met_dat[5] <- round(cor(stn_recon[time %in% genf$Jahr[genf$Jahr > 1900]], genf$Tagnr[genf$Jahr > 1900], use = "na.or.complete"),2)
    met_dat[6] <- round(mean((stn_recon[time %in% genf$Jahr[genf$Jahr < 1900]] - genf$Tagnr[genf$Jahr < 1900]), na.rm = T),2)
    met_dat[7] <- round(mean((stn_recon[time %in% genf$Jahr[genf$Jahr > 1900]] - genf$Tagnr[genf$Jahr > 1900]), na.rm = T),2)
    evallist_maesh13d_gve <- list(recon = stn_recon, obs = stn_obs, mets = met_dat)
  }

  if(grepl("mprua",names(ids[pv]))){
    liestalind <- c(which.min(abs(E - liestal_coord[1])),which.min(abs(N - liestal_coord[2])))
    stn_recon <- hcn[liestalind[1],liestalind[2],]
    stn_obs <- data.frame(reference_year = liestal$Jahr, doy = liestal$Tagnr)
    met_dat <- c()
    met_dat[1] <- round(cor(stn_recon[time %in% liestal$Jahr], liestal$Tagnr, use = "na.or.complete"),2)
    met_dat[2] <- round(mean((stn_recon[time %in% liestal$Jahr] - liestal$Tagnr)^2, na.rm = T)^0.5,2)
    met_dat[3] <- round(mean((stn_recon[time %in% liestal$Jahr] - liestal$Tagnr), na.rm = T),2)
    evallist_mprua65d_lti <- list(recon = stn_recon, obs = stn_obs, mets = met_dat)
  }

  if(grepl("mprua",names(ids[pv]))){
    ### swiss plateau mean
    nc <- nc_open(paste0(indir,"/",files_all[grepl("fldmean",files_all) & grepl("mprua",files_all)]))
    hcn <- ncvar_get(nc,"Mittelland")
    time <- lubridate::year(as.Date(ncvar_get(nc,"time"), origin = "1763-01-01"))
    stn_recon <- hcn
    stn_recon[1] <- NA
    stn_obs <- data.frame(reference_year = eurohist_dat$`Rutishauser, Kirschenbluete`$year,
                          doy = eurohist_dat$`Rutishauser, Kirschenbluete`$doy)
    met_dat <- c()
    met_dat[1] <- round(cor(stn_recon[time %in% stn_obs$reference_year], stn_obs$doy, use = "na.or.complete"),2)
    met_dat[2] <- round(mean((stn_recon[time %in% stn_obs$reference_year] - stn_obs$doy)^2, na.rm = T)^0.5,2)
    met_dat[3] <- round(mean((stn_recon[time %in% stn_obs$reference_year] - stn_obs$doy), na.rm = T),2)
    evallist_mprua65d_swp <- list(recon = stn_recon, obs = stn_obs, mets = met_dat)
  }
}


col1 = "#115f9a"; col2 = "#76c68f"; col3 = "#c9e52f"
lty <- 1; adj = 0.01
lwd1 <- 1.2; lwd2 <- 1.2; cex = 0.3; sigma = 3
alphaval <- 200

### load evaluation data #############3
for(pv in c("mprua65d","mfags13d","maesh13d")){
  (load(paste0("data/03_calibrations/",annealing_version,"/",pv,"/eval_mods_test_",pv,".RData")))
  assign(paste0(pv,"_test"), eval_mods_test)
  (load(paste0("data/03_calibrations/",annealing_version,"/",pv,"/eval_mods_train_",pv,".RData")))
  assign(paste0(pv,"_train"), eval_mods_train)
}

seedlen <- dim(mprua65d_test)[4]
mods <- colnames(eval_mods_test)
nc <- seedlen * length(mods)
cols <- rep(brewer.pal(n = 8, "YlGnBu"), each = seedlen)
nams <- c("maesh13d" = "Horse chestnut ","mprua65d" = "Cherry flowering","mfags13d" = "Beech leaf unfolding")
lwd = 0.4

rr <- list(c(0,25,5),c(-15,15,5),c(-0.3,1,0.2),c(-3,1,0.5))

metrics <- c("RMSE","Meanbias","Pcorr","NSE")

png(paste0(plotdir,"Figure_Phenoeval_obs_",Sys.Date(),".png"), width = 1200, height = 1200, res = 300, pointsize = 5)
par(mfcol = c(3,1), oma = c(3,1,3,1), cex = 1.2)
par(mar = c(0,4,0,0))
days <- seq(as.Date("1900-01-01"),as.Date("1900-07-01"), by = "day")
ind <- substr(days,9,10) %in% c("01","15");ind2 <- substr(days,9,10) %in% c("01")

plot(1763:2020, evallist_maesh13d_gve$recon, type = "l", ylab = "", xlab = "", yaxt = "n", xlim = c(1760, 2020), xaxt = "n",
     ylim = c(0,140), col = t_col(col1, alphaval), lwd = lwd1)
points(1763:2020, gaussian_wrap( evallist_maesh13d_gve$recon, sigma = 3), type = "l", col = col1, lwd = lwd1, lty = lty)
points(evallist_maesh13d_gve$obs$reference_year, evallist_maesh13d_gve$obs$doy, type = "l", col = t_col(col2, alphaval), lwd = lwd2, lty = lty, pch = 19, cex = cex)
points(evallist_maesh13d_gve$obs$reference_year, gaussian_wrap(evallist_maesh13d_gve$obs$doy, sigma = sigma), type = "l", col = col2, lwd = lwd2, lty = lty, pch = 19, cex = cex)
abline(h = which(ind2), lty = 2, col = "gray75", lwd = 0.5)
axis(side = 2, at = which(ind2) , substr(days,6,10)[ind2], las = 2)
mtext(side = 3, "Horse chestnut leaf unfolding", line  = -1.2, cex = 1.2, adj = adj)
mtext(side = 3, "b)", line  = 0.1, cex = 1.2, adj = adj)
legend("bottom", ncol = 2, col = c(col1, col2), legend = c("Reconstruction", "Observation"), lwd = lwd1, bty = "n")

plot(1763:2020, evallist_mprua65d_lti$recon, type = "l", ylab = "", xlab = "", yaxt = "n", xlim = c(1760, 2020), xaxt = "n",
     ylim = c(75,145), col = t_col(col1, alphaval), lwd = lwd1)
points(1763:2020, gaussian_wrap(evallist_mprua65d_lti$recon, sigma = sigma), type = "l", col = col1, lwd = lwd1 + 1, lty = lty)
points(evallist_mprua65d_lti$obs$reference_year, evallist_mprua65d_lti$obs$doy, type = "l", col = t_col(col2, alphaval),lwd = lwd2, lty = lty, pch = 19, cex = cex)
points(evallist_mprua65d_lti$obs$reference_year, gaussian_wrap(evallist_mprua65d_lti$obs$doy, sigma = sigma), type = "l", col = col2, lwd = lwd2, lty = lty, pch = 19, cex = cex)
abline(h = which(ind), lty = 2, col = "gray75", lwd = 0.5)
axis(side = 2, at = which(ind) , substr(days,6,10)[ind], las = 2)
mtext(side = 3, "Cherry tree flowering Liestal", line  = -1.2, cex = 1.2, adj = adj)

plot(1763:2020, evallist_mprua65d_swp$recon, type = "l", ylab = "", xlab = "", yaxt = "n", xlim = c(1760, 2020), xaxt = "n",
     ylim = c(75,145), col = t_col(col1, alphaval), lwd = lwd1)
points(1763:2020, gaussian_wrap(evallist_mprua65d_swp$recon, sigma = 3), type = "l", col = col1, lwd = lwd1 + 1, lty = lty)
points(eurohist_dat$`Rutishauser, Kirschenbluete`$year, eurohist_dat$`Rutishauser, Kirschenbluete`$doy, type = "l", col = t_col(col2, alphaval),lwd = lwd2, lty = lty, pch = 19, cex = cex)
points(eurohist_dat$`Rutishauser, Kirschenbluete`$year, gaussian_wrap(eurohist_dat$`Rutishauser, Kirschenbluete`$doy, sigma = sigma), type = "l", col = col2,lwd = lwd2, lty = lty, pch = 19, cex = cex)
abline(h = which(ind), lty = 2, col = "gray75", lwd = 0.5)
axis(side = 2, at = which(ind) , substr(days,6,10)[ind], las = 2)
mtext(side = 3, "Cherry tree flowering Swiss plateau", line  = -1.2, cex = 1.2, adj = adj)
axis(side = 1, at = seq(1760,2020,20) , seq(1760,2020,20), las = 1)

dev.off()

cols = c("#115f9a", "#76c68f", "#c9e52f", "#2a6395", "#8d5cb0", "#f7a71d", "#9d443b", "#3b8d6c")

### plot all the metrics in different plots
for(val in 1:4){
  png(paste0(plotdir,"Figure_Phenoeval_boxplots_allpvs_",metrics[val],"_",Sys.Date(),".png"),
      width = 1600, height = 1200, res = 300, pointsize = 5)
  par(mfrow = c(3,2), mar = c(0,0,0,0), oma = c(3,5,3,2), cex = 1.2) # c(3,1,3,1)

  for(pv in c("maesh13d","mprua65d","mfags13d")){
    test <- aperm(get(paste0(pv,"_test"))[,,val,], perm = c(1,3,2))
    train <- aperm(get(paste0(pv,"_train"))[,,val,], perm = c(1,3,2))
    ## change the dimensions later back!!!!!
    df <- data.frame(matrix(train,nrow =  nrow(train), ncol = prod(dim(train)[2:3]), byrow = F))
    boxplot(df[,3:ncol(df)], ylim = rr[[val]][1:2], xaxt = "n", col = cols[3:ncol(df)], outcex = 0.5, yaxt = "n")
    axis(side = 2, at = seq(rr[[val]][1],rr[[val]][2],rr[[val]][3]), las = 2)
    abline( h = seq(rr[[val]][1],rr[[val]][2],rr[[val]][3]), lty = 3, col = "gray75", lwd = lwd)
    if(pv == "mfags13d") axis(side = 1, at = c(1,seq(3,nc-3,3)), mods)
    if(pv == "maesh13d") {
      mtext(side = 3, text = "Training data", cex = 1.2)
      mtext(side = 3, text = "a)", adj = 0.01, cex = 1.2)}
    mtext(side = 2, text = nams[pv], line = 2.5, cex = 1.2)
    df <- data.frame(matrix(test,nrow = nrow(test), ncol =  prod(dim(train)[2:3]), byrow = F))
    boxplot(df[,3:ncol(df)], ylim = rr[[val]][1:2], xaxt = "n", col = cols[3:ncol(df)], yaxt = "n", outcex = 0.5)
    if(pv == "mfags13d") axis(side = 1, at = c(1,seq(3,nc-3,3)), mods)
    if(pv == "maesh13d") mtext(side = 3, text = "Testing data", cex = 1.2)
    abline( h = seq(rr[[val]][1],rr[[val]][2],rr[[val]][3]), lty = 3, col = "gray75", lwd = lwd)
  }
  dev.off()
}

