### HELPFER FUNCTIONS FOR COLD AND WARM SPELL CALCULATIONS #####

### parallelize quantile calculation voer grid
calc_quants_grid <- function(temp, date, qth){

  temp_list <- split(temp,seq(nrow(temp)))
  registerDoParallel(cores=floor(detectCores()*0.4))

  out <- foreach(ll = 1:length(temp_list),  .combine = "rbind") %dopar% {
    get_temp_quantiles(temp_rec = list(rec = temp_list[[ll]], dates = date),
                       qargs = list(qth=qth, n = 5, min_base_fraction=0.1,baseperiod=c(1871,1900), NAmaxQ=100))$base[,1]
  }
  return(out)
}

### quantile calculation over stations
calc_quants_station <- function(temp, date, qth){

  out <- get_temp_quantiles(list(rec = temp, dates = date),
                       qargs = list(qth=qth, n = 5, min_base_fraction=0.1,baseperiod=c(1871,1900), NAmaxQ=100))$base[,1]
  return(out)
}

### quantile function
get_temp_quantiles <- function (temp_rec,qargs)
{
  bootstrap_range <- as.Date(c(paste0(qargs$baseperiod,c("-01-01","-12-31"))))
  dates_base <-  seq(bootstrap_range[1], bootstrap_range[2], by = "day")

  temp <- temp_rec$rec
  base_val_nl <- which(as.character(temp_rec$dates) %in% as.character(dates_base) & substr(temp_rec$dates,6,10)!="02-29")
  temp_base <- temp[base_val_nl]

  outbase <- ClimIndVis:::zhang_running_qtile(temp_base, dates_base,qargs$qth, bootstrap_range = bootstrap_range,
                   n = qargs$n, min_fraction = qargs$min_base_fraction, NAmaxQ = qargs$NAmaxQ)
  return(list(base=outbase))

}


# get date factor from time
get_date_factors_obj_special<-function(t1d){
  factors<-list()
  factors$years<-factor(format(t1d,"%Y"))
  factors$yearmons<-factor(format(t1d,"%Y-%m"))
  jdays<-format(t1d,"%j")
  jdays<-replace_jday_29feb_new(jdays, t1d)
  factors$jdays=factor(jdays)
  return(factors)
}

get_date_factors_obj<-function(t1d){
  factors<-list()
  factors$years<-factor(format(t1d,"%Y"))
  factors$yearmons<-factor(format(t1d,"%Y-%m"))
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


replace_jday_29feb_special <- function(jdays){
  indices <- which(jdays == 366)
  if (length(indices) > 0){
    idays <- which(as.numeric(jdays) > 59)
    jdays[idays] <- pad2(as.numeric(jdays[idays]) - 1 , width=3)
  }
  return(jdays)
}


replace_jday_29feb_new <- function(jdays, t1d){
  if (any(substr(t1d,6,10) == "02-29")){
    iday <- which(substr(t1d,6,10) == "02-29")
    jdays[iday] <- jdays[iday - 1]
  }
  return(jdays)
}
