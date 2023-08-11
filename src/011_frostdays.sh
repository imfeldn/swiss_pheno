#!/bin/bash

temp_rec="../swiss_recon/temp/"

n=0
maxjobs=15

## frost days based on tmean
for file in $temp_rec/CH_temp_*.nc
 	
 	do 
 	
 	echo $file
	 . 011_frostdays_help.sh ${file} ${dir2} &
   
  if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished 
        echo $n wait
  fi
 
done

wait

echo "merge files"
cdo mergetime ${dir2}/fd/fd_yearly_* ${dir2}/fd/FD_yearly_1763-2020.nc
cdo ncrcat ${dir2}/fd/fd_yearly_* ${dir2}/fd/FD_yearly_1763-2020.nc
rm ${dir2}/fd/fd_yearly*
cdo ncrcat ${dir2}/fd/fd_monthly_* ${dir2}/fd/FD_monthly_1763-2020.nc
rm ${dir2}/fd/fd_monthly*
cdo mergetime ${dir2}/lfd/lfd_yearly_* ${dir2}/lfd/LFD_yearly_1763-2020.nc
rm ${dir2}/lfd/lfd_yearly_*
cdo mergetime ${dir2}/ffd/ffd_yearly_* ${dir2}/ffd/FFD_yearly_1763-2020.nc
rm ${dir2}/ffd/ffd_yearly_*

## calculate seasonal sums from monthly sums
cdo seassum ${dir2}/fd/FD_monthly_1763-2020.nc ${dir2}/fd/FD_season_1763-2020.nc

