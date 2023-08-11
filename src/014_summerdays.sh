#!/bin/bash

temp_rec="../swiss_recon/temp/"

# loop over years
for file in $temp_rec/CH_temp_*.nc
 	
 	do 
 	
 	# get years based on different data sets	
  if [[ $file == *"EnKF"* ]]; then year=$( echo "$file" | cut -c53-56); else year=$( echo "$file" | cut -c54-57); fi

 	echo $year	
        cdo yearsum -gec,25 $file ${dir2}/summerdays/summerdays_yearly_${year}.nc
  
        ## rename to make it consistent for later merging
        if [[ $file == *"TabsD"* ]]; then ncrename -v TabsD,temp -O ${dir2}/summerdays/summerdays_yearly_${year}.nc ${dir2}/summerdays/summerdays_yearly_${year}.nc; fi 

 
done

echo "merge files"

ncrcat ${dir2}/summerdays/summerdays_yearly_*.nc ${dir2}/summerdays/SUMMERDAYS_yearly_1763-2020.nc
cdo setctomiss,-999.99 ${dir2}/summerdays/SUMMERDAYS_yearly_1763-2020.nc ${dir2}/summerdays/SUMMERDAYS_yearly_1763-2020c.nc
mv ${dir2}/summerdays/SUMMERDAYS_yearly_1763-2020c.nc ${dir2}/summerdays/SUMMERDAYS_yearly_1763-2020.nc


