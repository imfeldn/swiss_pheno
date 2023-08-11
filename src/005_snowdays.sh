#!/bin/bash

temp_rec="../swiss_recon/temp/"
prec_rec="../swiss_recon/precip/"

start=1763
end=2020
years=($(seq ${start:0:4} ${end:0:4}))
len=${#years[@]}

n=0
maxjobs=10

# loop over years
for (( i = 0; i < ${len}; i++ ));do
 	
 	year=${years[${i}]}
 	 . 005_snowhelp.sh ${temp_rec} ${prec_rec} ${year} ${dir2} &
 
 if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        wait # wait until all have finished 
        echo $n wait
    fi
 
done

wait

echo "merge files"

cdo mergetime ${dir2}/snowdays/snow_monthly_* ${dir2}/snowdays/SD_monthly_${start}-${end}.nc
cdo mergetime ${dir2}/snowdays/snow_yearly_* ${dir2}/snowdays/SD_yearly_${start}-${end}.nc
cdo mergetime ${dir2}/snowdays/snow_season_* ${dir2}/snowdays/SNOWDAYS_season_${start}-${end}.nc

cdo mergetime ${dir2}/lsd/lsd_yearly_* ${dir2}/lsd/LSD_yearly_${start}-${end}.nc
cdo mergetime ${dir2}/fsd/fsd_yearly_* ${dir2}/fsd/FSD_yearly_${start}-${end}.nc

