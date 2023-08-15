###### SELECTION OF SERIES BY CLASS ##############

#### Packages
library(ncdf4)
library(stats)
library(geosphere)

#### directories
indir <- "wear/swiss_indices/"

### read station list
station_list <- read.csv(paste0(indir,"/2_data/pheno_net/data/ch.meteoschweiz.messnetz-phaenologie_de.csv"),header = T, sep = ",", allowEscapes = T)
row.names(station_list) <- as.character(station_list$Abk.)

### read all data
pheno <- read.table(paste0(indir,"/2_data/pheno_net/phaeno_previous.csv"),header=T, sep =";")

### read classification
classification <- read.table(paste0(indir,"/2_data/pheno_net/dataD4.csv"),header=T, sep =",")

### read qc
qc <- read.csv(paste0(indir,"/2_data/pheno_net/PhenoClass/QCflags_minimum_Final_Juni2017.csv"),
               header=T, sep =",")

### select stations with very valuable stations
#mprua65d	Cherry tree - full flowering 662
#maesh13d	Horse chestnut - leaf unfolding 601
#mfags13d	European beech - leaf unfolding 606

phenovals <- c("mfags13d", "mprua65d", "maesh13d","mcora65d","mcora13d","manen65d","mtaro65d","mtusf65d","mcarp65d","mlard13d")
param_id <- c(606,662,601,666,664,619,617,654,656,653,658,634)
names(param_id) <- phenovals

classes <- "123"

###### create files based on classes ############
for(pv in phenovals){

  vval_stns <- as.character(classification$station[which(classification$OLD_SHORT_TX == pv & classification$Rank %in% classes)])

  pheno_sel <- pheno[pheno$param_id == param_id[pv], ]
  vval_stns <- vval_stns[vval_stns %in% unique(pheno_sel$nat_abbr)]

  #length(unique(pheno$nat_abbr))
  pheno_sel <- pheno_sel[as.character(pheno_sel$nat_abbr) %in% vval_stns,]

  station_list_sel <- station_list[as.character(station_list$Abk.) %in% vval_stns,1:10]
  station_list_sel <- station_list_sel[vval_stns,] # reorder station list according to vvals

  # save as compressed RDS
  saveRDS(
    station_list_sel,
    file = paste0("data/02_pheno_net/stations_",pv,"_class",classes,".rds"),
    compress = "xz"
    )

  assign(paste0("station_list_",pv), station_list_sel)
  stnidx <- match(pheno_sel$nat_abbr, station_list_sel$Abk.)

  lat <- round(station_list_sel$Breitengrad[stnidx],5)
  lon <- round(station_list_sel$Laengengrad[stnidx],5)
  N <- station_list_sel$KoordinatenN[stnidx]
  E <- station_list_sel$KoordinatenE[stnidx]
  alt <- station_list_sel$Stationshoehe[stnidx]
  pheno_meta <- cbind(pheno_sel, lon, lat, E, N, alt)


  ### add qc flags
  qc_pv <- qc[qc$OLD_SHORT_TX == pv & qc$FinalFlagCountUnique == 1,c("station","year","doy","value","FinalFlagCountUnique")]

  for(ii in 1:nrow(qc_pv)){
    ind <- which(pheno_meta$nat_abbr == qc_pv[ii,"station"] & pheno_meta$reference_year == qc_pv[ii,"year"])
    if(length(ind) > 0){
      print(pheno_meta[ind,])
      pheno_meta[ind,"doy"] <- NA
    }
  }

  # save as compressed RDS
  saveRDS(
    pheno_meta,
    file = paste0("data/02_pheno_net/pheno_meta_",pv,"_class",classes,".rds"),
    compress = "xz"
    )

}
