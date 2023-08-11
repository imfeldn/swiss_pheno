#!/bin/bash


lpyears=(1764 1768 1772 1776 1780 1784 1788 1792 1796 1804 1808 1812 1816 1820 1824 1828 1832 1836 1840 1844 1848 1852 1856 1860 1864 1868 1872 1876 1880 1884 1888 1892 1896 1904 1908 1912 1916 1920 1924 1928 1932 1936 1940 1944 1948 1952 1956 1960 1964 1968 1972 1976 1980 1984 1988 1992 1996 2000 2004 2008 2012 2016 2020)

thmax=30
thmin=5

if [[ $file == *"EnKF"* ]]; then 
  
  year=$( echo "$file" | cut -c53-56); 
  # set values below/above threshold in file to thresholds
  cdo -expr,'temp = ((temp > 5)) ? temp : 5' $file ${dir2}/gdd/th5_${year}.nc
  cdo -expr,'temp = ((temp < 30)) ? temp : 30' ${dir2}/gdd/th5_${year}.nc ${dir2}/gdd/th530_${year}.nc

  else 
  year=$( echo "$file" | cut -c54-57)
  # set values below/above threshold in file to thresholds
  cdo -expr,'TabsD = ((TabsD > 5)) ? TabsD : 5' $file ${dir2}/gdd/th5_${year}.nc
  cdo -expr,'TabsD = ((TabsD < 30)) ? TabsD : 30' ${dir2}/gdd/th5_${year}.nc ${dir2}/gdd/th530_${year}.nc
  ncrename -v TabsD,temp -O ${dir2}/gdd/th530_${year}.nc ${dir2}/gdd/th530_${year}.nc

fi

echo ${year}
echo ${file}

if [[ " ${lpyears[*]} " =~ " ${year} " ]]; then doyfile="../doy_file_lp.nc"; else doyfile="../doy_file.nc"; fi
  	
# substract tbase and calculate cumulative sums
cdo timcumsum -subc,$thmin ${dir2}/gdd/th530_${year}.nc ${dir2}/gdd/gdd_${year}.nc
 		
# get days when a threshold is reached
cdo subc,1 -yearmin -add -setctomiss,0 -gec,1000 ${dir2}/gdd/gdd_${year}.nc $doyfile ${dir2}/gdd/gdd1000_${year}.nc
cdo subc,1 -yearmin -add -setctomiss,0 -gec,200 ${dir2}/gdd/gdd_${year}.nc $doyfile ${dir2}/gdd/gdd200_${year}.nc
cdo subc,1 -yearmin -add -setctomiss,0 -gec,1400 ${dir2}/gdd/gdd_${year}.nc $doyfile ${dir2}/gdd/gdd1400_${year}.nc
 	
rm ${dir2}/gdd/th5_${year}.nc
rm ${dir2}/gdd/th530_${year}.nc


