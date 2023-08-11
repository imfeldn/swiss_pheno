##### PLOT FOR APPENDIX PHENO EVALUATION BOXPLOTS ######
rm(list=ls())

### packages
library(RColorBrewer)

### select annleaing verision
annealing_version <- "sa_run_2023-08-09_class123"

for(pv in c("mprua65d","mfags13d","maesh13d")){
  (load(paste0("data/03_calibrations/",annealing_version,"/",pv,"/eval_mods_test_",pv,".RData")))
  assign(paste0(pv,"_test"), eval_mods_test)

  (load(paste0("data/03_calibrations/",annealing_version,"/",pv,"/eval_mods_train_",pv,".RData")))
  assign(paste0(pv,"_train"), eval_mods_train)
}

print("plot eval")

seedlen <- dim(mprua65d_test)[4]
mods <- colnames(eval_mods_test)
nc <- seedlen * ncol(eval_mods_test)
cols <- rep(brewer.pal(n = 8, "YlGnBu"), each = seedlen)
nams <- c("maesh13d" = "Horse chestnut ","mprua65d" = "Cherry flowering","mfags13d" = "Beech leaf unfolding")

metrics <- c("RMSE","Mean bias","Pearson correlation","NSE")
rr <- list(c(0,25,5),c(-20,20,5),c(-0.3,1,0.2),c(-3,1,0.5))
hmat <- matrix(1:6, nrow = 2, byrow = F); lwd <- 0.6

png(paste0("manuscript/supplementary/Figure_Appendix_Phenoeval_",Sys.Date(),".png"), width = 2000, height = 2300, res = 300, pointsize = 6)
layout(rbind(hmat, hmat + 6, hmat + 12))
par(oma = c(1,3,3,1))

for(pv in c("maesh13d","mprua65d","mfags13d")){

  for(val in 2:4){

    par(mar = c(0,2,2.5,1))
    test <- aperm(get(paste0(pv,"_test"))[,,val,], perm = c(1,3,2))
    train <- aperm(get(paste0(pv,"_train"))[,,val,], perm = c(1,3,2))
    df <- data.frame(matrix(train,nrow =  nrow(train), ncol = prod(dim(train)[2:3]), byrow = F))
    boxplot(df[,3:ncol(df)], ylim = rr[[val]][1:2], xaxt = "n", col = cols[3:ncol(df)], outcex = 0.5, yaxt = "n")
    axis(side = 2, at = seq(rr[[val]][1],rr[[val]][2],rr[[val]][3]), las = 2)
    abline( h = seq(rr[[val]][1],rr[[val]][2],rr[[val]][3]), lty = 3, col = "gray75", lwd = lwd)
    if(pv == "maesh13d") {mtext(side = 3, text = metrics[val], cex = 1, line = 2.5)}
    mtext(side = 3, text = nams[pv], cex = 1)
    if(val %in% 2) mtext(side = 2, text = "Training", cex = 1, line = 3)

    par(mar = c(2.5,2,0,1))
    df <- data.frame(matrix(test,nrow = nrow(test), ncol =  prod(dim(train)[2:3]), byrow = F))
    boxplot(df[,3:ncol(df)], ylim = rr[[val]][1:2], xaxt = "n", col = cols[3:ncol(df)], yaxt = "n", outcex = 0.5)
    if(val %in% 2)mtext(side = 2, text = "Testing", cex = 1, line = 3)
    abline( h = seq(rr[[val]][1],rr[[val]][2],rr[[val]][3]), lty = 3, col = "gray75", lwd = lwd)
    axis(side = 1, at = c(1,seq(3,nc-3,3)), mods)
    axis(side = 2, at = seq(rr[[val]][1],rr[[val]][2],rr[[val]][3]), las = 2)

  }
}
dev.off()
