#!/bin/bash


temp_rec="../swiss_recon/temp/"

n=0
maxjobs=6

for file in $temp_rec/CH_temp_*.nc
 	
 	do 
 	
  if [[ $file == *"EnKF"* ]]; then year=$( echo "$file" | cut -c53-56); else year=$( echo "$file" | cut -c54-57); fi

  echo $year
 	
  # get same missing values for all files
  cdo -f nc setmisstoc,-999.99 $file ${dir2}/tempmean/helpfile_${year}.nc
  
  # change the variable name
  if [[ $file == *"TabsD"* ]]
    then 
      ncrename -v TabsD,temp -O ${dir2}/tempmean/helpfile_${year}.nc ${dir2}/tempmean/helpfile_${year}.nc
  fi
 
  cdo monmean ${dir2}/tempmean/helpfile_${year}.nc ${dir2}/tempmean/tempmean_monthly_${year}.nc &
  cdo yearmean ${dir2}/tempmean/helpfile_${year}.nc ${dir2}/tempmean/tempmean_yearly_${year}.nc &
  
  if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished 
        echo $n wait
  fi

wait  	
 
done

rm ${dir2}/tempmean/helpfile_*
  
echo "merge files"

ncrcat ${dir2}/tempmean/tempmean_monthly_*.nc ${dir2}/tempmean/TEMPMEAN_monthly_1763-2020.nc
ncrcat ${dir2}/tempmean/tempmean_yearly_*.nc ${dir2}/tempmean/TEMPMEAN_yearly_1763-2020.nc
rm ${dir2}/tempmean/tempmean_monthly_*.nc
rm ${dir2}/tempmean/tempmean_yearly_*.nc


