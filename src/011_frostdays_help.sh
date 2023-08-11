#!/bin/bash

lpyears=(1764 1768 1772 1776 1780 1784 1788 1792 1796 1804 1808 1812 1816 1820 1824 1828 1832 1836 1840 1844 1848 1852 1856 1860 1864 1868 1872 1876 1880 1884 1888 1892 1896 1904 1908 1912 1916 1920 1924 1928 1932 1936 1940 1944 1948 1952 1956 1960 1964 1968 1972 1976 1980 1984 1988 1992 1996 2000 2004 2008 2012 2016 2020)

# frost days based on tmean 	
if [[ $file == *"EnKF"* ]]; then year=$( echo "$file" | cut -c53-56); else year=$( echo "$file" | cut -c54-57); fi

echo $year
newfile=${year}.nc
echo $newfile

if [[ " ${lpyears[*]} " =~ " ${year} " ]]; then doyfile="../doy_file_lp.nc"; else doyfile="../doy_file.nc"; fi
echo ${doyfile}
 	
cdo lec,0 ${file} ${dir2}/fd/le0_${newfile} 
## rename to make it consistent for later merging
if [[ $file == *"TabsD"* ]]; then ncrename -v TabsD,temp -O ${dir2}/fd/le0_${newfile} ${dir2}/fd/le0_${newfile}; fi 

## calculate monthly and yearly sums
cdo monsum ${dir2}/fd/le0_${newfile}  ${dir2}/fd/fd_monthly_${year}.nc
cdo yearsum ${dir2}/fd/le0_${newfile}  ${dir2}/fd/fd_yearly_${year}.nc

## last/first frost day
cdo setctomiss,0 ${dir2}/fd/le0_${newfile} ${dir2}/fd/le0_masked_${newfile}
cdo subc,1 -yearmax -selmonth,1/6 -add ${dir2}/fd/le0_masked_${newfile} $doyfile ${dir2}/lfd/lfd_yearly_${year}.nc
cdo subc,1 -yearmin -selmonth,7/12 -add ${dir2}/fd/le0_masked_${newfile} $doyfile ${dir2}/ffd/ffd_yearly_${year}.nc
	
rm ${dir2}/fd/le0_masked_${newfile}

