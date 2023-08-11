### HELPER FUNCTIONS TO CALCULATE CLIMATE INDICES #####

### growing season start
gsstart <- function(dat, date){

  if(is.list(date)) date <- date[[1]]

  f <- match.fun(">")
  min_length <- 6
  ind_mat <- aggregate(dat,list(year(date)), function(x){
    d <- f(x, 5)
    ## set to true if only one day in between is FALSE
    d[which(c(FALSE,d[1:length(d)-1]) & !d[1:length(d)] & c(d[2:length(d)],FALSE))] <- TRUE
    d2 <- Reduce(function(x, y) { # Reduce sucessively combines elements of a given vector
      z <- c(rep(FALSE, y), d[1:(length(d) - y)]) & x
      return(z)
    }, 1:(min_length-1), d)

    ## take the first one that gets 6 in a row and substract five.
    return(which(cumsum(d2) == 1)[1]) # in first version I substracted this five? return(which(cumsum(d2) == 1)[1] - 5)
  })
  colnames(ind_mat) <- c("year","annual")
  ind_mat <- as.matrix(ind_mat)
  return(ind_mat)
}

### growing season end
gsend <- function(dat, date){

  if(is.list(date)) date <- date[[1]]

  f <- match.fun("<")
  min_length <- 6
  ind_mat <- aggregate(dat,list(year(date)), function(x){
    d <- f(x, 5)
    ## set to true if only one day in between is FALSE
    d[which(c(FALSE,d[1:length(d)-1]) & !d[1:length(d)] & c(d[2:length(d)],FALSE))] <- TRUE
    d2 <- Reduce(function(x, y) { # Reduce sucessively combines elements of a given vector
      z <- c(rep(FALSE, y), d[1:(length(d) - y)]) & x
      return(z)
    }, 1:(min_length-1), d)
    return(which(cumsum(d2[182:length(d2)]) == 1)[1] + 182  -1)
  })

  colnames(ind_mat) <- c("year","annual")
  ind_mat <- as.matrix(ind_mat)
  return(ind_mat)
}

### growing season length
gsl <- function(dat, date){

  if(is.list(date)) date <- date[[1]]
  endind <- gsend(dat,date)[,"annual"]
  startind <- gsstart(dat,date)[,"annual"]
  gsl <- endind - startind
  ind_mat <- as.matrix(data.frame(unique(year(date)), gsl))
  colnames(ind_mat) <- c("year","annual")
  return(ind_mat)
}

# define gdd functions
tbase <- 5
tmax <- 30

gdd <- function(dat, date, th){

    if(is.list(date)){
      th <- date$th
      date <- date$date
    }

    dat[dat < tbase] <- tbase # no negative growing degree days
    dat[dat > tmax] <- tmax
    summed <- dat - tbase
    gdd <- aggregate(summed,list(year(date)), function(x) {
      y <- cumsum(x)
      if(any(is.na(x))){
        naind <- which(is.na(x))[1]-1
        if(naind != 0) if(y[naind] < th) {return(NA)}
      }
      if(!all(is.na(y))){
        if(any(round(y) >= th)){
          return(which.max(round(y) >= th)[1])
        } else{
          return(NA)
        }
      } else {
        return(NA)
      }
    })
    ind_mat <- as.matrix(gdd)
    colnames(ind_mat) <- c("year","annual")
    return(ind_mat)
}

### summerdays
summerdays <- function(dat, date){

  if(is.list(date)) date <- date[[1]]
  lastyr <- year(date[length(date)])

  seas <- get_seas_facts(date)
  ind <- (dat > 25) + 0
  ind_agg_seas <- aggregate(ind, list(seas), sum_Xpercna, na = 0.1)
  ind_agg_seas <- ind_agg_seas[as.numeric(substr(ind_agg_seas$Group.1,1,4)) <= lastyr,]
  ind_agg_yr <- aggregate(ind, list(year(date)), sum_Xpercna, na = 0.1)
  ind_mat <- matrix(ind_agg_seas$x, ncol = 4, nrow = length(unique(substr(ind_agg_seas$Group.1,1,4))), byrow = T)
  ind_mat <- cbind(ind_agg_yr$Group.1, ind_mat, ind_agg_yr$x)
  colnames(ind_mat) <- c("year", "DJF", "JJA","MAM","SON", "annual")
  return(ind_mat)
}

### frost days
fd <- function(dat, date){

  if(is.list(date)) date <- date[[1]]
  lastyr <- year(date[length(date)])

  seas <- get_seas_facts(date)
  ind <- (dat < 0) + 0
  ind_agg_seas <- aggregate(ind, list(seas), sum_Xpercna, na = 0.1)
  ind_agg_seas <- ind_agg_seas[as.numeric(substr(ind_agg_seas$Group.1,1,4)) <= lastyr,]
  ind_agg_yr <- aggregate(ind, list(year(date)), sum_Xpercna, na = 0.1)
  ind_mat <- matrix(ind_agg_seas$x, ncol = 4, nrow = length(unique(substr(ind_agg_seas$Group.1,1,4))), byrow = T)
  ind_mat <- cbind(ind_agg_yr$Group.1, ind_mat, ind_agg_yr$x)
  colnames(ind_mat) <- c("year", "DJF", "JJA","MAM","SON", "annual")
  return(ind_mat)
}

### first frost day
ffd <- function(dat, date){

  if(is.list(date)) date <- date[[1]]

  jdays <- get_date_factors_obj(date)$jdays
  ind <- ifelse((dat < 0),1,0) * as.numeric(jdays)
  ind_mat <- as.matrix(aggregate(ind, list(year(date)), function(x){
    halfyr <- x[round(365/2):length(x)]
    y <- halfyr[which(halfyr != 0)]
    z <- y[which.min(y)]
    if(any(which(is.na(halfyr)) < z)){z <- NA}
    if(length(z) == 0){z <- NA}
    return(z)
    }))
  colnames(ind_mat) <- c("year","annual")
  return(ind_mat)
}

### last frost day
lfd <- function(dat, date){

  if(is.list(date)) date <- date[[1]]

  jdays <- get_date_factors_obj(date)$jdays
  ind <- ifelse((dat < 0),1,0) * as.numeric(jdays)
  ind_mat <- as.matrix(aggregate(ind, list(year(date)), function(x){
    halfyr <- x[1:round(365/2)]
    if(all(is.na(halfyr))){
      return(NA)
      }else{
      y <- halfyr[which(halfyr != 0)]
      z <- y[which.max(y)]
      if(any(which(is.na(halfyr)) < z)){z <- NA}
      if(length(z) == 0){z <- NA}
      return(z)
    }
  }))
  colnames(ind_mat) <- c("year","annual")
  return(ind_mat)
}

### mean seasonal temperature
tempmean <- function(dat, date){

  if(is.list(date)) date <- date[[1]]
  lastyr <- year(date[length(date)])
  seanams <- c("DJF","JJA","MAM","SON")

  seas <- get_seas_facts(date)
  ind_agg_seas <- aggregate(dat, list(seas), mean_Xpercna, na = 0.1)
  ind_agg_yr <- aggregate(dat, list(year(date)), mean_Xpercna, na = 0.1)
  ind_mat <- matrix(NA, ncol = 6, nrow = nrow(ind_agg_yr), byrow = T)
  ind_mat[,1] <- ind_agg_yr$Group.1;ind_mat[,6] <- ind_agg_yr$x
  for(ss in 2:5){
    selyr <- which(substr(ind_agg_seas$Group.1,6,8) == seanams[ss-1] & as.numeric(substr(ind_agg_seas$Group.1,1,4)) <= lastyr)
    ind_mat[,ss]  <- ind_agg_seas$x[selyr]
    }
  colnames(ind_mat) <- c("year", seanams, "annual")
  return(ind_mat)
}


### create seasonal date factors
get_seas_facts <- function(dates){
  mds <- month(dates)
  seas_facts1 <- ifelse(mds %in% c(3:5),"MAM", ifelse(mds %in% c(6:8),"JJA",ifelse(mds %in% c(9:11),"SON","DJF")))
  seas_facts2 <- ifelse(mds == 12, year(dates) + 1, year(dates))
  seas_facts <- paste0(seas_facts2,"-",seas_facts1)
  return(as.factor(seas_facts))
}


### csdi
csdi <- function(temp, qdat, dates){

  if(is.list(dates)){dates <- dates$dates}
  lastyr <- year(dates[length(dates)])

  if(all(is.na(temp))){
    spelldays <- rep(NA,4)
    spellnumber<- rep(NA,4)
  }
  else{

    jdays <- get_date_factors_obj(dates)$jdays
    seas_facts <- get_seas_facts(dates)

    f <- match.fun("<")
    min_length <- 6
    d <- f(temp, qdat[jdays])
    ## set to true if only one day in between is FALSE
    d[which(c(FALSE,d[1:length(d)-1]) & !d[1:length(d)] & c(d[2:length(d)],FALSE))] <- TRUE

    d2 <- Reduce(function(x, y) { # Reduce sucessively combines elements of a given vector
      z <- c(rep(FALSE, y), d[1:(length(d) - y)]) & x
      return(z)
    }, 1:(min_length-1), d)

    periods <- Reduce(function(x, y) {
      return(c(d2[(y + 1):length(d2)], rep(FALSE, y)) | x)
    }, 1:(min_length-1), d2)

    spelldays  <- climdex.pcic:::tapply.fast(periods, seas_facts, sum_Xpercna, na = 0.1)

    times <- names(spelldays)
    yrs <- as.numeric(substr(times,1,4))
    spelldays <- spelldays[yrs <= lastyr]
    times <- times[yrs <= lastyr]
    yrs <- yrs[yrs <= lastyr]
  }

  ind_mat <- matrix(spelldays, ncol = 4, nrow = length(unique(yrs)), byrow = T)
  ind_mat <- cbind(unique(yrs), ind_mat)
  colnames(ind_mat) <- c("year", "DJF", "JJA","MAM","SON")

  return(ind_mat)
}


### wsdi
wsdi <- function(temp, qdat, dates){

  if(is.list(dates)){dates <- dates$dates}
  lastyr <- year(dates[length(dates)])

  if(all(is.na(temp))){
    spelldays <- rep(NA,4)
    spellnumber<- rep(NA,4)
  }
  else{
    jdays <- get_date_factors_obj(dates)$jdays
    seas_facts <- get_seas_facts(dates)

    f <- match.fun(">")
    min_length <- 6
    d <- f(temp, qdat[jdays])
    ## set to true if only one day in between is FALSE
    d[which(c(FALSE,d[1:length(d)-1]) & !d[1:length(d)] & c(d[2:length(d)],FALSE))] <- TRUE

    d2 <- Reduce(function(x, y) { # Reduce sucessively combines elements of a given vector
      z <- c(rep(FALSE, y), d[1:(length(d) - y)]) & x
      return(z)
    }, 1:(min_length-1), d)

    periods <- Reduce(function(x, y) {
      return(c(d2[(y + 1):length(d2)], rep(FALSE, y)) | x)
    }, 1:(min_length-1), d2)

    spelldays  <- climdex.pcic:::tapply.fast(periods, seas_facts, sum_Xpercna, na = 0.1)

    times <- names(spelldays)
    yrs <- as.numeric(substr(times,1,4))
    spelldays <- spelldays[yrs <= lastyr]
    times <- times[yrs <= lastyr]
    yrs <- yrs[yrs <= lastyr]
  }

  ind_mat <- matrix(spelldays, ncol = 4, nrow = length(unique(yrs)), byrow = T)
  ind_mat <- cbind(unique(yrs), ind_mat)
  colnames(ind_mat) <- c("year", "DJF", "JJA","MAM","SON")

  return(ind_mat)
}

monlist <- function(start_mon,end_mon){
  help=c(1:12,1:12)
  a<-which(help==start_mon)[1]
  b<-which(help == end_mon)
  x<-help[a:b[which(b>a)[1]]]
  return(x)
}


