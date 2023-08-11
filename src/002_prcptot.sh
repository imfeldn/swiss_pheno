#!/bin/bash

precip_rec="../swiss_recon/precip"

## loop over yearly files
for file in $precip_rec/CH_precip_*.nc; do 	
 	

 	if [[ $file == *"QMAP"* ]]; then year=$( echo "$file" | cut -c56-59); else year=$( echo "$file" | cut -c59-62); fi
 	echo ${year}
 	
 	cdo timsum ${file} ${dir2}/prcptot/prcptot_${year}.nc 
 	cdo monsum ${file}  ${dir2}/prcptot/prcptot_monthly_${year}.nc
 	
 	 # change the variable name
 	 if [[ $file == *"RhiresD"* ]]
  	  then 
  	    echo "test"
  	    ncrename -v RhiresD,precip -O ${dir2}/prcptot/prcptot_${year}.nc ${dir2}/prcptot/prcptot_${year}.nc
  	    ncrename -v RhiresD,precip -O ${dir2}/prcptot/prcptot_monthly_${year}.nc ${dir2}/prcptot/prcptot_monthly_${year}.nc
	  fi
 	
done

wait

echo "merge files"

cdo mergetime ${dir2}/prcptot/prcptot_* ${dir2}/prcptot/PRCPTOT_yearly_1763-2020.nc
cdo mergetime ${dir2}/prcptot/prcptot_monthly* ${dir2}/prcptot/PRCPTOT_monthly_1763-2020.nc


