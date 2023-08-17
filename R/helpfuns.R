#####  A COLLECTION OF HELPER FUNCTIONS #####

# extract leapyear from vector of years
leapyear <- function(x){
  return(ifelse((x%%4==0 & x%%100!=0) | x%%400==0,TRUE,FALSE))
}

# add padding
pad2 <- function(x, width=2){
  stringr::str_pad(x,width,pad="0")
}

# collapse
cllps <- function(x,pat=""){
  if(is.matrix(x)){
    sapply(1:nrow(x), function(z) paste0(x[z,],collapse = pat))
  }else{
    paste0(x,collapse = pat)
  }
}

# get transparent colors
t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}

# merge date columns from sef-files to date class
pdate <- function(dat){
	as.Date(paste0(dat[,1],"-",pad2(dat[,2]),"-",pad2(dat[,3])))
}

# calculate mean over values only when NA threshold is not reached
mean_Xpercna <- function(x, na=0.25){
  if(sum(is.na(x))/length(x) <= na){
    return(mean(x, na.rm=T))
  } else{return(NA)}}

# calculate sum over values only when NA threshold is not reached
sum_Xpercna <- function(x, na=0.25){
  if(sum(is.na(x))/length(x) <= na){
    return(sum(x, na.rm=T))
  } else{return(NA)}}

# change 29th february in doys
yday_without_29feb <- function(predperiod){
  leapyears <- leapyear(lubridate::year(predperiod))
  jind <- which(leapyears & substr(predperiod,6,10)>"02-28")
  jdays <- lubridate::yday(predperiod)
  if (length(jind) > 0) jdays[jind] <- jdays[jind]-1
  return(jdays)
}

# get window for quantile calculations
quant_window_fun <- function(dat,dates,doy){
  window <-8
  doysel <- (doy-window):(doy+window)
  doysel[which(doysel<1)] <- doysel[which(doysel<1)]+365
  doysel[which(doysel>365)] <- doysel[which(doysel>365)]-365
  doyind <- which(yday_without_29feb(dates) %in% doysel)
  q <- quantile(dat[doyind],probs=c(0.1,0.25,0.5,0.75,0.9), na.rm = T)
  return(q)
}

# remove seasonality of series based on two harmonics
rm_seas <- function(tempdata,dates) {

  if(!all(is.na(tempdata))){

    doys <- yday(dates)
    ndoy <- leapyear(year(dates))+365
    ndoy <- ndoy[1:length(dates)]
    xdoys <-doys/ndoy

    pcos <- cos(2*pi*xdoys)
    psin <- sin(2*pi*xdoys)
    pcos2 <- cos(4*pi*xdoys)
    psin2 <- sin(4*pi*xdoys)

    temp.jg.lm <- lm(tempdata ~ pcos + psin + pcos2 + psin2)
    temp.jg.cor <- as.numeric(tempdata[as.integer(names(temp.jg.lm$fitted.values))] - temp.jg.lm$fitted.values)
    tempdata[!is.na(tempdata)] <- temp.jg.cor
  }

  return(tempdata)
}


# get date factor from time
get_date_factors_obj<-function(t1d){
  factors<-list()
  factors$years<-factor(format(t1d,"%Y"))
  factors$yearmons<-factor(format(t1d,"%Y-%m"))
  factors$seas<-as.factor(ifelse(month(t1d) %in% c(3:5),"MAM", ifelse(month(t1d) %in% c(6:8),"JJA",ifelse(month(t1d) %in% c(9:11),"SON","DJF"))))
  jdays<-format(t1d,"%j")
  jdays<-replace_jday_29feb(jdays)
  factors$jdays=factor(jdays)
  return(factors)
}


replace_jday_29feb <- function(jdays){
  indices <- which(jdays == 366)
  if (length(indices) > 0)
    jdays[rep(indices, each = 366) + -365:0] <- pad2(c(1:59,59,60:365), width=3)
  return(jdays)
}


### evaluate reconstructions
eval_recon <- function(rec,obs){
  rmse <- mean((rec - obs)^2, na.rm = T)^0.5
  mae <- mean(abs(rec - obs),na.rm = T)
  pcorr <- cor(rec, obs, use = "pairwise.complete.obs")
  nse <- 1 - (sum((obs - rec)^2,na.rm = T)/sum((obs - mean(obs, na.rm = T))^2,na.rm = T))
  bias <- mean(rec - obs, na.rm = T)
  return(rbind(rmse,mae,pcorr,nse,bias))
}


