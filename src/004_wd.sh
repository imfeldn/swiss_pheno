#!/bin/bash

precip_rec="/scratch3/noemi/wear/swiss_recon/precip/"

n=0
maxjobs=10

# loop over yearly files
for file in $precip_rec/CH_precip_*.nc; do 
 	
  
  if [[ $file == *"QMAP"* ]]; then year=$( echo "$file" | cut -c57-60); else year=$( echo "$file" | cut -c60-63); fi
  echo ${year}
  
 	cdo gec,1 ${file} ${dir2}/wd/precip_gec01_${year}.nc 

  
  # change the variable name
  if [[ $file == *"RhiresD"* ]]
    then 
      ncrename -v RhiresD,precip -O ${dir2}/wd/precip_gec01_${year}.nc ${dir2}/wd/precip_gec01_${year}.nc
  fi
  
  cdo monsum ${dir2}/wd/precip_gec01_${year}.nc  ${dir2}/wd/wd_monthly_${year}.nc &
  cdo yearsum ${dir2}/wd/precip_gec01_${year}.nc  ${dir2}/wd/wd_yearly_${year}.nc &
   
  if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished 
        echo $n wait
    fi

done

wait

rm ${dir2}/wd/precip_gec01_*.nc 

echo "merge files"

ncrcat ${dir2}/wd/wd_yearly* ${dir2}/wd/WD_yearly_1763-2020.nc
ncrcat ${dir2}/wd/wd_monthly* ${dir2}/wd/WD_monthly_1763-2020.nc

cdo setctomiss,-999.99 ${dir2}/wd/WD_yearly_1763-2020.nc ${dir2}/wd/WD_yearly_1763-2020.nc
cdo setctomiss,-999.99 ${dir2}/wd/WD_monthly_1763-2020.nc ${dir2}/wd/WD_monthly_1763-2020.nc
cdo seassum ${dir2}/wd/WD_monthly_1763-2020.nc ${dir2}/wd/WD_season_1763-2020.nc


