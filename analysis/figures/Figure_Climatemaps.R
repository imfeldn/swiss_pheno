########  PLOT OF ALL CLIMATE INDICES #####
rm(list=ls())

### packages
library(ncdf4);library(lubridate)
library(viridis);library(RColorBrewer)

### functions
source("R/helpfuns.R")
source("R/image_funs.R")

### directories
plotdir <- "manuscript/figures/"

### indices to display
indices_spring <- list("gsstart"=list(var="temp",tdef="yearly"), "gdd" = list(var="temp",tdef="yearly",agg = "200"),
                       "lfd"=list(var="temp",tdef="yearly"),"fd" = list(var="temp",tdef="season"),
                       "wsdi" = list(var="days",tdef="season"),"csdi" = list(var="days",tdef="season"),
                       "snowdays" = list(var="temp",tdef="season"),"wd"=list(var="precip",tdef="season"))

nams <- names(indices_spring)
longnams <- c("Growing season start" ,"Growing degree days","Last frost day","Frost day","Warm spell index","Cold spell index","Snow days","Wet days")

### make a list of the seasonal indices
spring_list <- list()
tt <- 0

for(ii in names(indices_spring)){

  print(ii)
  files <- list.files(path = paste0("data/01_climate_indices/",ii,"/"), pattern = indices_spring[[ii]]$tdef, full.names = T)
  files <- files[!grepl("fldmean|regions|old|hwnumber|cwnumber",files)]
    tt <- tt + 1

    if(ii == "gdd"){
      file_sel <- grepl(indices_spring[[ii]]$agg,files)
      nc <- nc_open(files[file_sel])
      dat <- ncvar_get(nc,indices_spring[[ii]]$var)
      dat[dat < -999] <- NA
      time <- as.Date(ncvar_get(nc,"time"),origin = "1970-01-01")
      spring_list[[tt]]  <- list(dat= dat,time = time)
      nc_close(nc)
    } else{
      nc <- nc_open(files)
      dat <- ncvar_get(nc,indices_spring[[ii]]$var)
      dat[dat < -999] <- NA
      time <- as.Date(ncvar_get(nc,"time"),origin = "1970-01-01")
      spring_list[[tt]]  <- list(dat= dat,time = time)
      nc_close(nc)
    }
  }

names(spring_list) <- nams
rm(dat)

### calc climatology for each period
periods <- list(per1 = c(1781,1810), per2 = c(1811,1840), per3 = c(1841,1870) ,per_ind = c(1871,1900),
                 per5 = c(1901,1930), per6 = c(1931,1960), per7 = c(1961,1990), per8 = c(1991,2020))

dd <- dim(spring_list$gsstart$dat)[1:2]
dd2 <- dim(spring_list$wd$dat)[1:2]

na_th <- function(x){sum(is.na(x))/length(x) > 0.8}
calc <- F

if(calc){
  for (ind in 1:length(spring_list)){
    print(nams[ind])
    dx <- switch((ind == 8) + 1,dd,dd2)
    ## calculate climatology
    climmat <- array(sapply(periods, function(x){
      xind <- which(year(spring_list[[ind]]$time) %in% c(x[1]:x[2]))
      return(apply(spring_list[[ind]]$dat[,,xind],1:2,mean, na.rm = T))
    }), dim = c(dx,length(periods)))
    assign(paste0("clim_",nams[ind]), climmat)
    save(climmat, file = paste0("data/01_climate_indices/clims/clim_",nams[ind],".RData"))

    ## calculate t-test
    pmat <- array(sapply(periods, function(p){
      print(p)
      xind <- which(year(spring_list[[ind]]$time) %in% c(p[1]:p[2]))
      xind_pre <- which(year(spring_list[[ind]]$time) %in% c(periods$per_ind[1]:periods$per_ind[2]))
      xdat <- plyr::alply(spring_list[[ind]]$dat[,,xind],c(1,2))
      ydat <- plyr::alply(spring_list[[ind]]$dat[,,xind_pre],c(1,2))
      mat <- mapply(function(x,y){if (any(na_th(x),na_th(y))){return(NA)} else{t.test(x,y, na.action ="na.pass")$p.value} },x = xdat, y = ydat)
      return(array(data = mat, dim = dx))
    }), dim = c(dx,length(periods)))
    assign(paste0("pval_",nams[ind]), pmat)
    save(pmat, file = paste0("data/01_climate_indices/clims/pval_",nams[ind],".RData"))
  }

} else{

  for(ind in 1:length(spring_list)){
    load(paste0("data/01_climate_indices/clims/clim_",nams[ind],".RData"));assign(paste0("clim_",nams[ind]), climmat)
    load(paste0("data/01_climate_indices/clims/pval_",nams[ind],".RData"));assign(paste0("pval_",nams[ind]), pmat)
  }
}

rm(climmat,pmat)

####### plot all indices at once ####################
brks <- list("gsstart" = seq(40,175,10),"fd" = seq(0,95,5), "lfd" = seq(30,170,10), "csdi" = seq(0,15,1),"snowdays" = seq(0,30,3),"wd" = seq(30,62,4),
             "gdd" = seq(90,230,10), "wsdi" = seq(5,35,2.5))
cols <- list("gsstart" = colorRampPalette(rev(c("#6F5741","#B48222","#E1D075","#EBF3B1","#CCFF9C","#C9E959","#9ACD32","#4FBD33","#619E40","#326329")))(length(brks[["gsstart"]])-1),
             "gdd" = colorRampPalette(rev(c("#6F5741","#B48222","#E1D075","#EBF3B1","#CCFF9C","#C9E959","#9ACD32","#4FBD33","#619E40","#326329")))(length(brks[["gdd"]])-1),
             "fd" = colorRampPalette(c("#0E6134","#2BAC66","#ABCF63","#E8F69E","#D7E3EE","#B5CAFF","#8FB3FF","#7F97FF","#0570B0","#023858"))(length(brks[["fd"]])-1),
             "lfd" = colorRampPalette(c("#0E6134","#2BAC66","#ABCF63","#E8F69E","#D7E3EE","#B5CAFF","#8FB3FF","#7F97FF","#0570B0","#023858"))(length(brks[["lfd"]])-1),
             "csdi" = colorRampPalette(c("#0E6134","#2BAC66","#ABCF63","#E8F69E","#D7E3EE","#B5CAFF","#8FB3FF","#7F97FF","#0570B0","#023858"))(length(brks[["csdi"]])-1),
             "wsdi" = colorRampPalette(c("#FFFF8F","#FED679","#FBBA5B","#F99C4F","#E76636","#D62F27","#B50626","#6B0026","#2B0016"))(length(brks[["wsdi"]])-1),
             "wd" = colorRampPalette(c("#EDFAC2","#CDFFCD","#99F0B2","#53BD9F","#32A696","#3296B4","#0570B0","#05508C","#0A1F96","#2C0246","#6A2C5A"))(length(brks[["wd"]])-1),
             "snowdays" = colorRampPalette(c("#EDFAC2","#CDFFCD","#99F0B2","#53BD9F","#32A696","#3296B4","#0570B0","#05508C","#0A1F96","#2C0246","#6A2C5A"))(length(brks[["snowdays"]])-1)
)

indsel <- c(1,2,5,3,4,6,7,8)

marscale <- c(2,1,2,5);cx <- 1.5;adj <- 0.5;cx1 <- 1.5
png(paste0(plotdir,"/Figure_Climatemaps_",Sys.Date(),".png"), width = length(periods) * 700, height = length(indsel) * 550, res = 300, pointsize = 7)
layout(matrix(1:((length(periods) + 1)*length(indsel)), nrow = length(indsel), ncol = length(periods) + 1, byrow = T), widths = c(rep(1,length(periods)),0.26), heights = 1)
par(oma = c(1,2,2,1));it = 0
for(isel in indsel){
  par(mar = c(1,1,1,1));it = it + 1
  dat <- get(paste0("clim_",nams[isel]))
  dat[dat < brks[[nams[isel]]][1]] <- brks[[nams[isel]]][1]
  dat[dat > brks[[nams[isel]]][length(brks[[nams[isel]]])]] <- brks[[nams[isel]]][length(brks[[nams[isel]]])]
  for(i in 1:length(periods)) {
    image(1:nrow(dat),1:ncol(dat),dat[,,i], xaxt  = "n", yaxt = "n", bty = "n", breaks = brks[[nams[isel]]], col = cols[[nams[isel]]])
    if(isel == 1) mtext(side = 3, text = paste0(periods[[i]][1]," - ", periods[[i]][2]), adj = adj, cex  = cx1)
    if(i == 1) mtext(side = 2, longnams[isel], cex = cx1)
    if(i == 1) mtext(side = 3, paste0(letters[it],")"), cex = cx1, adj = 0.01)
  }
  image_scale(breaks=brks[[nams[isel]]],col=cols[[nams[isel]]], axis.pos=4, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
              useraxis = list(axis.pos=4,ats=brks[[nams[isel]]][seq(1,length(brks[[nams[isel]]]),2)],labs=brks[[nams[isel]]][seq(1,length(brks[[nams[isel]]]),2)], cex=cx,axis.pos2=1,line=3, adj=-0.07))
}
dev.off()


####### plot all indices at once ANOMALIES ####################
brks <- list("gsstart" = seq(-24,24,4),"fd" = seq(-10,10,2), lfd = seq(-25,25,5), csdi = seq(-10,10,2), wsdi = seq(-20,20,4),"gdd" = seq(-25,25,5),
             "snowdays" = seq(-5,5,1),"wd" = seq(-5,5,1))
cols <- list("gsstart" = colorRampPalette(rev(c("#6F5741","#B48222","white","#619E40","#326329")))(length(brks[["gsstart"]])-1),
             "gdd" =colorRampPalette(rev(c("#6F5741","#B48222","white","#619E40","#326329")))(length(brks[["gdd"]])-1),
             "fd" = brewer.pal(n = length(brks[["fd"]])-1, name = "RdBu"),
             "lfd" = brewer.pal(n = length(brks[["lfd"]])-1, name = "RdBu"),
             "csdi" = brewer.pal(n = length(brks[["csdi"]])-1, name = "RdBu"),
             "wsdi" = rev(brewer.pal(n = length(brks[["wsdi"]])-1, name = "RdBu")),
             "wd" = colorRampPalette(c("#6F5741","#B48222","white","#619E40","#326329"))(length(brks[["wd"]])-1),
             "snowdays" = colorRampPalette(c("#6F5741","#B48222","white","#619E40","#326329"))(length(brks[["snowdays"]])-1))

selper <- c(1:length(periods))
png(paste0(plotdir,"/Figure_Climatemaps_Anomalies_",Sys.Date(),".png"),  width = length(periods) * 700, height = length(indsel) * 550, res = 300, pointsize = 7)
layout(matrix(1:((length(periods) + 1)*length(indsel)), nrow = length(indsel), ncol = length(periods) + 1, byrow = T), widths = c(rep(1,length(periods)),0.26), heights = 1)
par(oma = c(1,2,2,1));it = 0
for(isel in c(1:8)){
  par(mar = c(1,1,1,1));it = it + 1
  dat <- sweep(get(paste0("clim_",nams[isel])),MARGIN = c(1,2),STATS = get(paste0("clim_",nams[isel]))[,,names(periods)=="per_ind"])
  dat[dat < brks[[nams[isel]]][1]] <- brks[[nams[isel]]][1]
  dat[dat > brks[[nams[isel]]][length(brks[[nams[isel]]])]] <- brks[[nams[isel]]][length(brks[[nams[isel]]])]
  for(i in selper) {
    image(1:nrow(dat),1:ncol(dat),dat[,,i], xaxt  = "n", yaxt = "n", bty = "n", breaks = brks[[nams[isel]]], col = cols[[nams[isel]]])
    if(isel == 1) mtext(side = 3, text = paste0(periods[[i]][1]," - ", periods[[i]][2]), adj = adj, cex  = cx1)
    if(i == 1) mtext(side = 2, longnams[isel], cex = cx1)
    if(i == 1) mtext(side = 3, paste0(letters[it],")"), cex = cx1, adj = 0.01)
  }
  image_scale(breaks=brks[[nams[isel]]],col=cols[[nams[isel]]], axis.pos=4, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
              useraxis = list(axis.pos=4,ats=brks[[nams[isel]]][seq(1,length(brks[[nams[isel]]]),2)],labs=brks[[nams[isel]]][seq(1,length(brks[[nams[isel]]]),2)], cex=cx,axis.pos2=1,line=3, adj=-0.07))
}
dev.off()

####### plot all indices at once significance ####################
brks_sig <- c(0,0.05,seq(0.1,1, 0.1))
cols_sig <- colorRampPalette(c("#225ea8", "#41b6c4", "#a1dab4", "#ffffcc", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026"))(length(brks_sig)-1)

png(paste0(plotdir,"/Figure_Climatemaps_Pvalues_",Sys.Date(),".png"),  width = length(periods) * 700, height = length(indsel) * 550, res = 300, pointsize = 7)
layout(matrix(1:((length(selper)+1) *length(indices_spring)), nrow = length(indices_spring), ncol = length(selper) + 1, byrow = T), widths = c(rep(1,length(selper)),0.27), heights = 1)
par(oma = c(1,2,2,1));it = 0
for(isel in c(1:length(indices_spring))){
  par(mar = c(1,1,1,1));it = it + 1
  dat <- get(paste0("pval_",nams[isel]))
  dat2 <- dat
  dat2[dat >= 0.05] <- NA
   for(i in selper) {
    image(1:nrow(dat),1:ncol(dat),dat[,,i], xaxt  = "n", yaxt = "n", bty = "n", breaks = brks_sig, col = cols_sig)
    #image(1:nrow(dat),1:ncol(dat),dat2[,,i], xaxt  = "n", yaxt = "n", bty = "n", col = "pink", add = T)
    if(isel == 1) mtext(side = 3, text = paste0(periods[[i]][1]," - ", periods[[i]][2]), adj = adj, cex  = cx1)
    if(i == 1) mtext(side = 2, longnams[isel], cex = cx1)
    if(i == 1) mtext(side = 3, paste0(letters[it],")"), cex = cx1, adj = 0.01)
  }
  image_scale(breaks=brks_sig,col=cols_sig, axis.pos=4, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
              useraxis = list(axis.pos=4,ats=brks_sig[seq(1,length(brks_sig),2)],labs=brks_sig[seq(1,length(brks_sig),2)], cex=cx,axis.pos2=1,line=3, adj=-0.07))
}
dev.off()
