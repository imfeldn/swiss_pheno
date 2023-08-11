##### MODEL EVALUATION #######
rm(list=ls())

### load packages
library(phenor);library(RColorBrewer)

### functions
source("R/helpfuns.R")

AICc_nna <- function (measured, predicted, k) {
  n <- length(measured)
  RSS <- sum((measured - predicted)^2, na.rm = T)
  AIC <- 2 * k + n * log(RSS/n)
  AICc <- AIC + (2 * k * (k + 1))/(n - k - 1)
  return(list(AIC = AIC, AICc = AICc))
}

### select subdir to evaluate
subDir <- "sa_run_2023-08-09_class123"

### select phenophases for calibration
# phenovals <- c("mfags13d","mprua65d","maesh13d",
#                "mtaro65d","manen65d","mmald65d",
#                "mpyrc65d","mcora65d","mcora13d",
#                "mtusf65d","mcarp65d","mlard13d")
phenovals <- c("maesh13d","mfags13d","mprua65d")

#pv <- "mprua65d"

### evaluate phenophases
for(pv in phenovals){

  print(pv)

  # load model parameters, testing and training data
  load(paste0("data/03_calibrations/",subDir,"/",pv,"/model_comparison_",pv,".RData"))
  load(paste0("data/03_calibrations/",subDir,"/",pv,"/stn_list_train_",pv,".RData"))
  load(paste0("data/03_calibrations/",subDir,"/",pv,"/stn_list_test_",pv,".RData"))

  seeds <- 1:nrow(mod_out$modelled$LIN$parameters)
  mods <- names(mod_out$modelled)

  # evaluation training
  eval_mods_train <- array(NA, dim=c(length(stn_list_train),length(mods),4,length(seeds)))
  dimnames(eval_mods_train) <- list(names(stn_list_train),mods, c("rmse","bias","pcorr","nse"), seeds)

  for (mod in 1:length(mod_out$modelled)){

    for(ss in 1:length(seeds)){

      t <- 0
      for(stn in 1:length(stn_list_train)){
        stnind <- (t + 1):(t + length(stn_list_train[[stn]]$transition_dates))
        t <- stnind[length(stnind)]
        obs <- stn_list_train[[stn]]$transition_dates
        pred <- mod_out$modelled[[mod]]$predicted_values[ss,]
        eval_mods_train[stn,mod,,ss] <- round(eval_recon(rec = pred[stnind], obs)[c("rmse","bias","pcorr","nse"),],3)
      }
    }
  }
  save(eval_mods_train, file = paste0("data/03_calibrations/",subDir,"/",pv,"/eval_mods_train_",pv,".RData"))

  # evaluation testing

  stn_list_test <-  stn_list_test[sapply(stn_list_test, function(x) length(x$year)) > 1]

  eval_mods_test <- array(NA, dim=c(length(stn_list_test),length(mods),4,length(seeds)))
  dimnames(eval_mods_test) <- list(names(stn_list_test),mods, c("rmse","bias","pcorr","nse"),seeds)

  for (mod in 1:length(mod_out$modelled)){

    #print(mod)
    for(ss in 1:length(seeds)){

      pred <- pr_predict(par = mod_out$modelled[[mod]]$parameters[ss,], data = stn_list_test, model = names(mod_out$modelled)[mod])
      pred[pred > 1000] <- NA
      t <- 0

      for(stn in 1:length(stn_list_test)){
        stnind <- (t + 1):(t + length(stn_list_test[[stn]]$transition_dates))
        t <- stnind[length(stnind)]
        obs <- stn_list_test[[stn]]$transition_dates
        eval_mods_test[stn,mod,,ss] <- round(eval_recon(rec = pred[stnind], obs)[c("rmse","bias","pcorr","nse"),],3)
        }
    }
  }
  save(eval_mods_test, file = paste0("data/03_calibrations/",subDir,"/",pv,"/eval_mods_test_",pv,".RData"))

  ### calculate AICc
  aicc_mat <- array(NA, dim = c(length(mods),length(seeds)));rownames(aicc_mat) <- mods

  for(mm in 1:length(mods)){
    k <- ncol(mod_out$modelled[[mm]]$parameters)
    aicc_mat[mm,] <- sapply(1:length(seeds), function(x) {AICc_nna(measured = mod_out$measured, predicted = mod_out$modelled[[mm]]$predicted_values[x,], k = k)$AICc})
  }

  save(aicc_mat, file = paste0("data/03_calibrations/",subDir,"/",pv,"/aicc_mods_",pv,".RData"))

  # overall model output
  # feval <- sapply(mod_out$modelled, function(x){
  #   dat <- x$predicted_values
  #   return(sapply(1:length(seeds), function(y) round(eval_recon(dat[y,],mod_out$measured),3)))
  # })
  # rownames(feval) <- rep(c("rmse","mae","corr","nse","bias"), length(seeds))
  # write.table(feval, file = paste0("data/calibrations/",subDir,"/",pv,"/meaneval_mods_train_",pv,".txt"))

  # plot the evaluation ###
  print("plot eval")

  nrtest <- nrow(eval_mods_test)
  nrtrain <- nrow(eval_mods_train)
  nc <- length(seeds) * length(mods)
  cols <- rep(brewer.pal(n = length(mods), "YlGnBu"), each = length(seeds))

  train <- aperm(eval_mods_train, perm = c(1,4,3,2))
  test <- aperm(eval_mods_test, perm = c(1,4,3,2))

  lwd = 0.5
  png(paste0("data/03_calibrations/",subDir,"/",pv,"/testtrain_boxplots_",pv,"_",Sys.Date(),".png"), width = 3200, height = 1300, res = 300, pointsize = 6)
  par(mfcol = c(2,4), mar = c(1,2,0,1), oma = c(2,2,2,2))
  rr <- c(floor(min(c(test[,,1,], train[,,1,]), na.rm = T)),ceiling(max(c(test[,,1,], train[,,1,]), na.rm = T)))
  if (rr[2] > 25) rr[2] <- 25
  df <- data.frame(matrix(train[,,1,],nrow = nrtrain, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  mtext(side = 3,line = 0.1, "RMSE")
  abline( h = seq(0,25,5), lty = 3, col = "gray75", lwd = lwd)
  df <- data.frame(matrix(test[,,1,],nrow = nrtest, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  axis(side = 1, at = seq(2,nc,length(seeds)), mods)
  abline( h = seq(0,25,5), lty = 3, col = "gray75", lwd = lwd)

  rr <- c(floor(min(c(test[,,2,], train[,,2,]), na.rm = T)),ceiling(max(c(test[,,2,], train[,,2,]), na.rm = T)))
  df <- data.frame(matrix(train[,,2,],nrow = nrtrain, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  mtext(side = 3,line = 0.1, "Mean bias")
  abline( h = seq(-25,25,5), lty = 3, col = "gray75", lwd = lwd)
  df <- data.frame(matrix(test[,,2,],nrow = nrtest, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  axis(side = 1, at = seq(2,nc,length(seeds)), mods)
  abline( h = seq(-25,25,5), lty = 3, col = "gray75", lwd = lwd)

  rr <- c(floor(min(c(test[,,3,], train[,,3,]), na.rm = T)),ceiling(max(c(test[,,3,], train[,,3,]), na.rm = T)))
  df <- data.frame(matrix(train[,,3,],nrow = nrtrain, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  mtext(side = 3,line = 0.1, "Pearson correlation")
  abline( h = seq(-1,1,0.2), lty = 3, col = "gray75", lwd = lwd)
  df <- data.frame(matrix(test[,,3,],nrow = nrtest, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  axis(side = 1, at = seq(2,nc,length(seeds)), mods)
  abline( h = seq(-1,1,0.2), lty = 3, col = "gray75", lwd = lwd)

  rr <- c(floor(min(c(test[,,4,], train[,,4,]), na.rm = T)),ceiling(max(c(test[,,4,], train[,,4,]), na.rm = T)))
  if (rr[1] < -3) rr[1] <- -3
  df <- data.frame(matrix(train[,,4,],nrow = nrtrain, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  mtext(side = 3,line = 0.1, "NSE")
  abline( h = seq(-4,1,0.5), lty = 3, col = "gray75", lwd = lwd)
  df <- data.frame(matrix(test[,,4,],nrow = nrtest, ncol = nc, byrow = F))
  boxplot(df, ylim = rr, xaxt = "n", col = cols)
  axis(side = 1, at = seq(2,nc,length(seeds)), mods)
  abline( h = seq(-8,1,0.5), lty = 3, col = "gray75", lwd = lwd)
  dev.off()

}
