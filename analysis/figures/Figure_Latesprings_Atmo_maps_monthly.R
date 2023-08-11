##### PLOT LATE SPRING INDICES #######
rm(list=ls())

### packages
library(ncdf4);library(lubridate);library(RColorBrewer)

### helpfer functions
source("R/image_funs.R")
source("R/plot_winds.R")
#source("../../ncode/R/randomfuns.R")

### directores
plotdir <- "manuscript/figures/"

latesprgs <- c(1785,1837,1853)

### changes this path still ? #####
(load("../swiss_indices/5_output/4_pheno_recon/Springindex_phenorecon_2023-06-13.RData"))

early <- df_springs_yr$years[order(df_springs_yr$ranks)][1:10]
late <- rev(df_springs_yr$years[order(df_springs_yr$ranks)])[2:10]

latesprgs <- late#c(1785,1837,1853)

################ modera late springs ################
files <- list.files("../swiss_indices/2_data/modera/",pattern = "europe", full.names = T)
files <- files[grepl("geopot|slp|temp|totprec|_u_|_v_", files) & !grepl("season", files)]

## read files
nc <- nc_open(files[1])
thelp <- paste0(rep(1421:2009, each = 12),"-",pad2(1:12,2), "-15")
nc_time <- as.Date(thelp[ncvar_get(nc, "time") + 1])
years <- unique(year(nc_time))
lats <- rev(ncvar_get(nc,"latitude"))
lons <- ncvar_get(nc,"longitude")

for(ff in 1:length(files)){

  mylist <- list()

  nc <- nc_open(files[ff])
  (varname <- names(nc$var)[1])
  mylist[["lats"]] <- rev(ncvar_get(nc,"latitude"))
  mylist[["lons"]] <- ncvar_get(nc,"longitude")
  mylist[["nc_time"]] <- as.Date(thelp[ncvar_get(nc, "time") + 1])
  mylist[["years"]] <- unique(year(nc_time))
  mat <- ncvar_get(nc,varname)
  mat <- mat[,ncol(mat):1,]
  mylist[["mat"]] <- array(mat, dim = c(dim(mat)[1:2],12,dim(mat)[3]/12))
  mylist[["climmat"]] <- apply(mylist[["mat"]][,,,years %in% c(1871:1900)],1:3,mean)
  mylist[["anommat"]] <- sweep(mylist[["mat"]],MARGIN = 1:3, STATS = mylist[["climmat"]])
  assign(paste0(varname,"_list"), mylist)
  nc_close(nc)
}

rm(mylist)

brksp <- seq(-200,200,25)
colp <- rev(colorRampPalette(brewer.pal(n =5, name = "RdBu"))(length(brksp)-1))


###### the three latest springs in one plot ##########
marscale = c(1,0.5,1,3);cx <- 1
selsprings <- c(1785,1837,1853)

png(paste0(plotdir,"Figure_Latesprings_Atmos_all_",Sys.Date(), ".png"), width = 3300, height = 2200, res = 300, pointsize = 10)
layout(matrix(c(1:15), nrow = 3, ncol = 5, byrow = T), heights = c(1), widths = c(1,1,1,1,0.22))
par(oma = c(2,2,2,1))

for(ll in 1:length(selsprings)){
  par(mar = c(1,0.5,1,0.5))
  ## modera upper-level
  for(mm in 1:4){
    mat = geopotential_height_list$anommat[,,mm + 1,years == selsprings[ll]]
    mat[mat < brksp[1]] = brksp[1]
    mat[mat > brksp[length(brksp)]] = brksp[length(brksp)]
    xaxt <- ifelse(mm %in% c(2,4),"n", "n")
    yaxt <- ifelse(mm %in% c(1,2),"n", "n")
    image(geopotential_height_list$lons,geopotential_height_list$lats,mat, col = colp, breaks = brksp, xaxt = xaxt, yaxt = yaxt)
    maps::map("world", add = T, col = "gray55")
    maps::map("world",region = "Switzerland", add = T, col = "red")
    contour(geopotential_height_list$lons,geopotential_height_list$lats,geopotential_height_list$mat[,,mm + 1,years == selsprings[ll]],
            xaxt = "n", yaxt = "n", add = T, lwd = 0.6, drawlabels = F)
    seqlon = seq(1,length(eastward_wind_list$lons),4)
    seqlat = seq(1,length(eastward_wind_list$lats),4)

    quiver(u=eastward_wind_list$mat[seqlon,seqlat,mm + 1,years == selsprings[ll]],v=northward_wind_list$mat[seqlon,seqlat,mm + 1,years == selsprings[ll]],
           x.val= eastward_wind_list$lons[seqlon],y.val= eastward_wind_list$lats[seqlat],scale=6, length=1, col="black",lwd=0.5)

    if(mm == 1) mtext(side = 3, adj = 0.01, paste0(letters[ll],")"), line = 0.2)
    mtext(side = 3, paste0(c("February","March","April","May")[mm]," ",selsprings[ll]), line = 0.1)
  }
  image_scale(breaks=brksp,col=colp, axis.pos=4, scale_lab="",  add.axis=F,mar=marscale,key.extend = c(F,F),
              useraxis = list(axis.pos=4,ats=brksp[seq(1,length(brksp),2)],labs=brksp[seq(1,length(brksp),2)], cex=1,axis.pos2=4,line=2, adj=0.01))
}

dev.off()

