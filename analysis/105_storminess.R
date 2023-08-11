#### SLP STANDARD DEVIATION AFTER APPLYING A LANCZOS PASS-BAND FILTER  #####

### packages
library(lubridate);library(XLConnect)

springs <- c(1770, 1785, 1816, 1837, 1853, 1932)

### Lanczos filter
Lf <- t(read.table("analysis/lanczos.txt"))

lanczos_filter <- function(x){
  if(any(is.na(x))){
    return(NA)
  } else{
    sdval <- sd(filter(x, Lf, method = "convolution"), na.rm=TRUE)
    return(sdval)
  }
}

### load data
load("../swiss_arm/code/data_input/analogue_stationdata.RData")

### read inventory
inventory <- readWorksheetFromFile("../swiss_arm/code/data_input/arm_stations_inventory_withnbcn.xlsx",
                                   sheet=1, colTypes="character",startRow=2)

### calculate bandpass filtered sd climatology of the stations for 1961- 1990
pind <- which(grepl("_p",colnames(TOT)) & !grepl("_prob",colnames(TOT)))
TOTp <- TOT[,c(1,pind)]

yrs <- unique(year(TOT$date))
sd_yearly_MAM <- matrix(NA, ncol = ncol(TOTp) - 1, nrow = length(yrs))
rownames(sd_yearly_MAM) <- 1763:2020
colnames(sd_yearly_MAM) <- names(TOTp)[2:ncol(TOTp)]

## select spring days before and after
seldays <- c(46:166)

for(yr in 1:length(yrs)){
  tind <- yday(TOTp$date) %in% seldays & year(TOTp$date) == yrs[yr]
  sd_yearly_MAM[yr,] <- apply(TOTp[tind,2:ncol(TOTp)],2, lanczos_filter)
}

climyrs <- 1961:1990#:2020
climind <- yrs %in% climyrs
sd_clim_MAM <- colMeans(sd_yearly_MAM[climind,], na.rm = T)
sd_anom_MAM <- sweep(sd_yearly_MAM, MARGIN = c(2), STATS = sd_clim_MAM, FUN = "-")

save(file = paste0("analysis/storminess_MAM_lanzcos_ref19611990_",Sys.Date(),".RData"), sd_anom_MAM)
