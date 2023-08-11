#!/bin/bash

lpyears=(1764 1768 1772 1776 1780 1784 1788 1792 1796 1804 1808 1812 1816 1820 1824 1828 1832 1836 1840 1844 1848 1852 1856 1860 1864 1868 1872 1876 1880 1884 1888 1892 1896 1904 1908 1912 1916 1920 1924 1928 1932 1936 1940 1944 1948 1952 1956 1960 1964 1968 1972 1976 1980 1984 1988 1992 1996 2000 2004 2008 2012 2016 2020)

if [[ $file == *"EnKF"* ]]; then year=$( echo "$file" | cut -c53-56); else year=$( echo "$file" | cut -c54-57); fi
#year=$( echo "$file" | cut -c53-56)
echo $year

if [[ " ${lpyears[*]} " =~ " ${year} " ]]
 	then 
 		doyfile="../doy_file_lp.nc"
 		lastday=366
 	else 
 		doyfile="../doy_file.nc"
 		lastday=365
fi
echo ${doyfile}
 
# start 
cdo gtc,5 $file ${dir2}/gsl/mask_$year.nc
cdo consecsum ${dir2}/gsl/mask_$year.nc ${dir2}/gsl/periods_$year.nc
cdo eqc,6 ${dir2}/gsl/periods_$year.nc ${dir2}/gsl/eq6periods_$year.nc
cdo setctomiss,0 ${dir2}/gsl/eq6periods_$year.nc ${dir2}/gsl/help.nc 
cdo subc,1 -yearmin -selmonth,1/6 -add ${dir2}/gsl/help.nc $doyfile ${dir2}gsl/gsstart_${year}.nc
rm ${dir2}/gsl/mask_$year.nc ${dir2}/gsl/periods_$year.nc ${dir2}/gsl/eq6periods_$year.nc ${dir2}/gsl/help.nc
 	
# end  
cdo ltc,5 $file ${dir2}/gsl/mask_$year.nc
cdo consecsum ${dir2}/gsl/mask_$year.nc ${dir2}/gsl/periods_$year.nc
cdo eqc,6 ${dir2}/gsl/periods_$year.nc ${dir2}/gsl/eq6periods_$year.nc
cdo setctomiss,0 ${dir2}/gsl/eq6periods_$year.nc ${dir2}/gsl/help.nc 
cdo subc,1 -yearmin -selmonth,7/12 -add ${dir2}/gsl/help.nc $doyfile ${dir2}gsl/gsend_${year}.nc
rm ${dir2}/gsl/mask_$year.nc ${dir2}/gsl/periods_$year.nc ${dir2}/gsl/eq6periods_$year.nc ${dir2}/gsl/help.nc

# length
cdo sub ${dir2}gsl/gsend_${year}.nc ${dir2}gsl/gsstart_${year}.nc ${dir2}gsl/gsl_${year}.nc

# change the variable name
if [[ $file == *"TabsD"* ]]
  then 
  ncrename -v TabsD,temp -O ${dir2}gsl/gsend_${year}.nc ${dir2}gsl/gsend_${year}.nc
  ncrename -v TabsD,temp -O ${dir2}gsl/gsstart_${year}.nc ${dir2}gsl/gsstart_${year}.nc
  ncrename -v TabsD,temp -O ${dir2}gsl/gsl_${year}.nc ${dir2}gsl/gsl_${year}.nc
fi


