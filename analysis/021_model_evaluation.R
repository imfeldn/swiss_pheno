##### MODEL EVALUATION #######
rm(list=ls())

### load packages
library(phenor);library(RColorBrewer)

### functions
source("R/helpfuns.R")

AICc_nna <- function (measured, predicted, k) {
  message("only use if you compare same measurements")
  nnaind <- which(!is.na(measured))
  n <- length(measured[nnaind])
  RSS <- sum((measured[nnaind] - predicted[nnaind])^2)
  AIC <- 2 * k + n * log(RSS/n)
  AICc <- AIC + (2 * k * (k + 1))/(n - k - 1)
  return(list(AIC = AIC, AICc = AICc))
}

### select subdir to evaluate
classes <- "123"
subDir <- paste0("data/03_calibrations/sa_run_2023-08-15_class",classes)

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
  mod_out <- readRDS(paste0(subDir,"/",pv,"/model_comparison_",pv,".rds"))
  stn_list_train <- readRDS(paste0(subDir,"/",pv,"/stn_list_train_",pv,".rds"))
  stn_list_test <- readRDS(paste0(subDir,"/",pv,"/stn_list_test_",pv,".rds"))
  train_sites <- unique(stn_list_train$site)
  test_sites <- unique(stn_list_test$site)

  seeds <- 1:nrow(mod_out$modelled$LIN$parameters)
  mods <- names(mod_out$modelled)

  # evaluation training
  eval_mods_train <- array(NA, dim=c(length(stn_list_train),length(mods),4,length(seeds)))
  dimnames(eval_mods_train) <- list(names(stn_list_train),mods, c("rmse","bias","pcorr","nse"), seeds)

  for (mod in 1:length(mod_out$modelled)){
    for(ss in 1:length(seeds)){
      for(stn in 1:length(stn_list_train)){
        stnind <- stn_list_train$site == train_sites[stn]
        obs <- stn_list_train$transition_dates[stnind]
        pred <- mod_out$modelled[[mod]]$predicted_values[ss,stnind]
        eval_mods_train[stn,mod,,ss] <- round(eval_recon(rec = pred, obs)[c("rmse","bias","pcorr","nse"),],3)
      }
    }
  }

  saveRDS(
    eval_mods_train,
    file = file.path(subDir,pv,paste0("model_evaluation_train_",pv,".rds")),
    compress = "xz"
  )
  # evaluation testing
  eval_mods_test <- array(NA, dim=c(length(stn_list_test),length(mods),4,length(seeds)))
  dimnames(eval_mods_test) <- list(names(stn_list_test),mods, c("rmse","bias","pcorr","nse"),seeds)

  for (mod in 1:length(mod_out$modelled)){

    for(ss in 1:length(seeds)){

      pred <- pr_predict(par = mod_out$modelled[[mod]]$parameters[ss,], data = stn_list_test, model = names(mod_out$modelled)[mod])
      pred[pred > 1000] <- NA ## set to NA

      for(stn in 1:length(stn_list_test)){
        stnind <- stn_list_test$site == test_sites[stn]
        obs <- stn_list_test$transition_dates[stnind]
        eval_mods_test[stn,mod,,ss] <- round(eval_recon(rec = pred[stnind], obs)[c("rmse","bias","pcorr","nse"),],3)
        }
    }
  }

  saveRDS(
    eval_mods_test,
    file = file.path(subDir,pv,paste0("model_evaluation_test_",pv,".rds")),
    compress = "xz"
  )

  ### calculate AICc
  aicc_mat <- array(NA, dim = c(length(mods),length(seeds)));rownames(aicc_mat) <- mods

  for(mm in 1:length(mods)){
    k <- ncol(mod_out$modelled[[mm]]$parameters)
    aicc_mat[mm,] <- sapply(1:length(seeds), function(x) {AICc_nna(measured = mod_out$measured, predicted = mod_out$modelled[[mm]]$predicted_values[x,], k = k)$AICc})
  }

  saveRDS(
    aicc_mat,
    file = file.path(subDir,pv,paste0("model_evaluation_test_",pv,".rds")),
    compress = "xz"
  )

  # plot the evaluation ###
  print("plot eval")

  nrtest <- nrow(eval_mods_test)
  nrtrain <- nrow(eval_mods_train)
  nc <- length(seeds) * length(mods)
  cols <- rep(brewer.pal(n = length(mods), "YlGnBu"), each = length(seeds))

  train <- aperm(eval_mods_train, perm = c(1,4,3,2))
  test <- aperm(eval_mods_test, perm = c(1,4,3,2))

  lwd = 0.5
  png(paste0("data/03_calibrations/",subDir,"/",pv,"/testtrain_boxplots_",pv,"_",Sys.Date(),".png"), width = 3800, height = 1300, res = 300, pointsize = 6)
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
