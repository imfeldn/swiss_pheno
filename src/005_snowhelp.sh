#!/bin/bash


temp_rec=${temp_rec}
precip_rec=${precip_rec}
echo ${year}
  
# select appropriate files  
temp_file=${temp_rec}/CH_temp_*_${year}*.nc
prec_file=${prec_rec}/CH_precip_*_${year}*.nc

lpyears=(1764 1768 1772 1776 1780 1784 1788 1792 1796 1804 1808 1812 1816 1820 1824 1828 1832 1836 1840 1844 1848 1852 1856 1860 1864 1868 1872 1876 1880 1884 1888 1892 1896 1904 1908 1912 1916 1920 1924 1928 1932 1936 1940 1944 1948 1952 1956 1960 1964 1968 1972 1976 1980 1984 1988 1992 1996 2000 2004 2008 2012 2016 2020)

if [[ " ${lpyears[*]} " =~ " ${year} " ]]; then doyfile="../doy_file_lp.nc"; else doyfile="../doy_file.nc"; fi
echo ${doyfile}

echo $temp_file  	

if [[ ${year} > 1960 ]]
  then 
  echo "change var name"
  ncrename -v TabsD,temp -O ${temp_file} ${dir2}/snowdays/helpfile_${year}.nc
  cdo ltc,2 ${dir2}/snowdays/helpfile_${year}.nc ${dir2}/snowdays/lt2_temp_${year}.nc
  rm ${dir2}/snowdays/helpfile_${year}.nc
  
  else 
  cdo ltc,2 ${temp_file} ${dir2}/snowdays/lt2_temp_${year}.nc
fi

cdo selindexbox,1,370,6,245 -gec,1 ${prec_file} ${dir2}/snowdays/gec1_precip_${year}.nc
cdo eqc,2 -add ${dir2}/snowdays/lt2_temp_${year}.nc ${dir2}/snowdays/gec1_precip_${year}.nc ${dir2}/snowdays/sd_${year}.nc
 	
cdo monsum ${dir2}/snowdays/sd_${year}.nc  ${dir2}/snowdays/snow_monthly_${year}.nc 
cdo yearsum ${dir2}/snowdays/sd_${year}.nc  ${dir2}/snowdays/snow_yearly_${year}.nc 

cdo setctomiss,0 ${dir2}/snowdays/sd_${year}.nc ${dir2}/snowdays/sd0_masked_${year}.nc
cdo subc,1 -yearmax -selmonth,1/6 -add ${dir2}/snowdays/sd0_masked_${year}.nc $doyfile ${dir2}/lsd/lsd_yearly_${year}.nc
cdo subc,1 -yearmin -selmonth,7/12 -add ${dir2}/snowdays/sd0_masked_${year}.nc $doyfile ${dir2}/fsd/fsd_yearly_${year}.nc

rm ${dir2}/snowdays/lt2_temp_${year}.nc;rm ${dir2}/snowdays/gec1_precip_${year}.nc;rm ${dir2}/snowdays/sd_${year}.nc
 
