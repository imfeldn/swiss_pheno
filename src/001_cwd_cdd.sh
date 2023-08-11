#!/bin/bash

#set -o errexit
#set -o pipefail

precip_rec="../swiss_recon/precip/"

for file in $precip_rec/CH_precip_*.nc
 	
 	do 
 	 
        if [[ $file == *"QMAP"* ]]; then year=$( echo "$file" | cut -c57-60); else year=$( echo "$file" | cut -c60-63); fi
        echo ${year}

 	cdo eca_cwd $file ${dir2}/cwd/cwd_${year}.nc
 	cdo eca_cdd $file ${dir2}/cdd/cdd_${year}.nc
	
done

echo "merge files"

cdo mergetime ${dir2}/cwd/cwd_* ${dir2}/cwd/CWD_1763-2020.nc
rm ${dir2}/cwd/cwd_*
cdo mergetime ${dir2}/cdd/cdd_* ${dir2}/cdd/CDD_1763-2020.nc
rm ${dir2}/cdd/cdd_*

