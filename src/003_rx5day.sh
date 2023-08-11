#!/bin/bash


precip_rec="../swiss_recon/precip/"


for file in $precip_rec/CH_precip_R*.nc; do 

 	if [[ $file == *"QMAP"* ]]; then year=$( echo "$file" | cut -c57-60); else year=$( echo "$file" | cut -c60-63); fi
 	echo ${year}
 	
 	cdo timmax -ydrunsum,5 $file ${dir2}/rx5day/rx5day_${year}.nc 
  
  	 # change the variable name
 	 if [[ $file == *"RhiresD"* ]]
  	  then 
  	    echo "test"
  	    ncrename -v RhiresD,precip -O ${dir2}/rx5day/rx5day_${year}.nc ${dir2}/rx5day/rx5day_${year}.nc
    fi


done

echo "merge files"

wait

ncrcat ${dir2}/rx5day/rx5day_* ${dir2}/rx5day/RX5DAY_yearly_1763-2020.nc
cdo setctomiss,-999.99 ${dir2}/rx5day/RX5DAY_yearly_1763-2020.nc ${dir2}/rx5day/RX5DAY_yearly_1763-2020c.nc
mv ${dir2}/rx5day/RX5DAY_yearly_1763-2020c.nc ${dir2}/rx5day/RX5DAY_yearly_1763-2020.nc


