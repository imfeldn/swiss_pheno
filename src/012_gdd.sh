#!/bin/bash

# check the following link out
#https://code.mpimet.mpg.de/boards/2/topics/11131

temp_rec="../swiss_recon/temp/"

n=0
maxjobs=20

# loop over years
for file in $temp_rec/CH_temp_*.nc
 	
 	do 
 	
	 . 012_gdd_help.sh ${file} ${dir2} &
   
  if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished 
        echo $n wait
  fi

wait  	
 
done

echo "merge files"


cdo ncrcat gdd/gdd1000_* gdd/GDD1000_yearly_1763-2020.nc
rm gdd/gdd1000_*
cdo ncrcat gdd/gdd200_* gdd/GDD200_yearly_1763-2020.nc
rm gdd/gdd200_*



