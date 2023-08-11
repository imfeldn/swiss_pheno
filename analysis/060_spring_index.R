##### CREATE SPRING INDEX FOR ENTIRE PERIOD #####

### packages
library(ncdf4);library(lubridate);library(smoother)

#### functions ??
source("/scratch3/noemi/wear/swiss_arm/code/013_helpfun.R")
source("R/helpfuns.R")

# PCA of incex???
# https://towardsdatascience.com/a-one-stop-shop-for-principal-component-analysis-5582fb7e0a9c
doys <- substr(seq(as.Date("1990-01-01"),as.Date("1990-12-31"), by = "day"),6,10)[1:360]

### read topo in order to only use low elevation variables
nc <- nc_open("../swiss_arm/code/data_input/topo.swiss02_ch01r.swisscors.nc")
topo <- ncvar_get(nc, "height")
topo2000 <- which(topo < 1600)

### read Regulas spring index
sprg_mch <- read.table("../swiss_indices/2_data/spring_index/pheno_phenospringindex_2022_G.txt", skip = 14, header = T, sep = "\t")

### read the pheno reconstructions
path <- "data/04_pheno_recon/"
files <- list.files(path = path, pattern = c("CH_"))#maesh|mprua|mfags|mcora|mtaro|mcarp|manen|mlard|mtusf
files <- files[!grepl("swissplateau|index|fldmean", files)]
pvs <- gsub("CH_M1s*.|CH_M1*.|CH_PTTs*.|CH_PTT*.","",files)
pvs <- substr(pvs,1,8)

ff <- 1

for(ff in 1:length(files)){
  nc <- nc_open(paste0(path, files[ff]))
  pheno <- ncvar_get(nc, "DOY")
  pheno[,,1] <- NA
  pheno_time <- year(as.Date(ncvar_get(nc,"time"), origin = "1763-01-01"))
  pheno_ts <- sapply(1:length(pheno_time), function(x){
    y <- pheno[,,x]
    return(mean(y[topo2000], na.rm = T))
  })
  assign(paste0(pvs[ff],"_ts"),pheno_ts)
  nc_close(nc)
  rm(pheno,pheno_ts)
}

mpheno_ts <- (mfags13d_ts + maesh13d_ts + mprua65d_ts)/3
plot(pheno_time, mpheno_ts, type = "l", lwd = 2)
points(sprg_mch$year, sprg_mch$index + 115, type = "l", col = "forestgreen", lwd = 2)
mtext(side = 3, "Spring index comparison")

cor1 = round(cor.test(sprg_mch$index[sprg_mch$year %in% pheno_time], mpheno_ts[pheno_time %in% sprg_mch$year])$estimate,2)

df_mod <- data.frame(pred = sprg_mch$index[sprg_mch$year %in% pheno_time],
                 mfags13d = mfags13d_ts[pheno_time %in% sprg_mch$year],
                 mprua65d = mprua65d_ts[pheno_time %in% sprg_mch$year],
                 maesh13d = maesh13d_ts[pheno_time %in% sprg_mch$year],
                 mcora13d = mcora13d_ts[pheno_time %in% sprg_mch$year],
                 mcora65d = mcora65d_ts[pheno_time %in% sprg_mch$year],
                 mcarp65d = mcarp65d_ts[pheno_time %in% sprg_mch$year],
                 mtaro65d = mtaro65d_ts[pheno_time %in% sprg_mch$year],
                 manen65d = manen65d_ts[pheno_time %in% sprg_mch$year],
                 mlard13d = mlard13d_ts[pheno_time %in% sprg_mch$year],
                 mtusf65d = mtusf65d_ts[pheno_time %in% sprg_mch$year]
                 )

df_hist <- data.frame(mfags13d = mfags13d_ts,
                      mprua65d = mprua65d_ts,
                      maesh13d = maesh13d_ts,
                      mcora13d = mcora13d_ts,
                      mcora65d = mcora65d_ts,
                      mcarp65d = mcarp65d_ts,
                      mtaro65d = mtaro65d_ts,
                      manen65d = manen65d_ts,
                      mlard13d = mlard13d_ts,
                      mtusf65d = mtusf65d_ts
                      )

#define intercept-only model
intercept_only <- lm(pred ~ 1, data=df_mod)
summary(intercept_only)

#define model with all predictors
all <- lm(pred ~ ., data=df_mod)
summary(all)

#perform backward stepwise regression
stepfun <- step(intercept_only, direction='both', scope=formula(all), trace=3)
summary(stepfun)

threespecies <- lm(pred ~ mfags13d + maesh13d + mprua65d, data=df_mod)
summary(threespecies )

spring_pred1 <- predict(all,df_hist)
spring_pred2 <- predict(stepfun,df_hist)
spring_pred3 <- predict(threespecies,df_hist)

View(data.frame(1763:2020,round(spring_pred1,2),round(spring_pred2,2),round(spring_pred3,2)))

round(eval_recon(rec = spring_pred2[pheno_time %in% sprg_mch$year], obs = sprg_mch$index[sprg_mch$year %in% pheno_time]),3)
round(eval_recon(rec = spring_pred1[pheno_time %in% sprg_mch$year], obs = sprg_mch$index[sprg_mch$year %in% pheno_time]),3)
round(eval_recon(rec = spring_pred3[pheno_time %in% sprg_mch$year], obs = sprg_mch$index[sprg_mch$year %in% pheno_time]),3)


plot(pheno_time, round(spring_pred1,2), type = "l", col = "forestgreen", lwd = 2.5, xlab = "", ylab = "", yaxt = "n", ylim = c(-15,20), xlim = c(1950,2020))
points(pheno_time, round(spring_pred2,2), type = "l", col = "cornflowerblue", lwd = 1)
points(pheno_time, round(spring_pred3,2), type = "l", col = "magenta", lwd = 1)
points(sprg_mch$year, sprg_mch$index, type = "l", col = "gold2", lwd = 2)
axis(side = 2, at = seq(-20,20,5), labels = seq(-20,20,5), col = "forestgreen", las = 2)


png(paste0(plotdir,"Spring_index_recon_stepwise_", Sys.Date(),".png"), width = 800, height = 300)
plot(pheno_time, mpheno_ts, type = "l", lwd = 2.5, ylim = c(80, 150), xlab = "", ylab = "", yaxt = "n", col = "gray35")
axis(side = 2, at = seq(1,250, 15), labels = doys[seq(1,250, 15)], las = 2)
abline(h = seq(1,250, 15), lty = 3, col = "gray75")
par(new = T)
plot(pheno_time, spring_pred2, type = "l", col = "forestgreen", lwd = 2.5, xlab = "", ylab = "", yaxt = "n", ylim = c(-15,40))
points(sprg_mch$year, sprg_mch$index, type = "l", col = "gold2", lwd = 2)
axis(side = 4, at = seq(-12,12,4), labels = seq(-12,12,4), col = "forestgreen", las = 2)
abline(h = 0, lty = 3, col = "forestgreen")
mtext(side = 3, "Spring Index Reconstruction 1763 - 2020")
legend("bottom",ncol = 2, bty = "n", lwd = 2,
       col = c("black","forestgreen","gold", "cornflowerblue"),
       legend = c("mean of phenos","lin.reg of spring index", "MCH spring index (no unit)", "index swiss plateau"))
dev.off()

#spring_class <- ifelse(vals < -6.58, 1,ifelse(vals < -3.46, 2,ifelse(vals < 3.46, 3,ifelse(vals < 6.58, 4, 5))))
df_springs_yr <- data.frame(years = 1763:2020,
                            index = round(spring_pred2,4),
                            mean_pheno = round(mpheno_ts,2),
                            date  = c(doys[round(mpheno_ts)]))
View(df_springs_yr[order(df_springs_yr$index),])

save(df_springs_yr, file =  paste0("data/04_pheno_recon/springindex_phenorecon_", Sys.Date(),".RData"))
#save(sprg_mch, file =  paste0("data/04_pheno_recon/Springindex_MCH_", Sys.Date(),".RData"))

