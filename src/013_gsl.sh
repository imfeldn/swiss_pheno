#!/bin/bash

#### this is not adjusted at all!! #####

temp_rec="../swiss_recon/temp/"

n=0
maxjobs=19

for file in $temp_rec/CH_temp_*.nc
 	
 	do 
 	
  . 013_gsl_help.sh $file &
  
  if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished 
        echo $n wait
  fi

wait  	
 
done

echo "merge files"

ncrcat ${dir2}/gsl/gsl_*.nc ${dir2}/gsl/GSL_yearly_1763-2020.nc
rm ${dir2}/gsl/gsl_*.nc
ncrcat ${dir2}/gsl/gsstart_*.nc ${dir2}/gsl/GSSTART_yearly_1763-2020.nc
rm ${dir2}/gsl/gsstart_*.nc
ncrcat ${dir2}/gsl/gsend_*.nc ${dir2}/gsl/GSEND_yearly_1763-2020.nc
rm ${dir2}/gsl/gsend_*.nc

