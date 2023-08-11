#!/bin/bash

set -o errexit
set -o pipefail

dir2="data/climate_indices/"


echo "calc frost indices"
. 011_frostdays.sh ${dir2} 
echo "calc gdd"
. 012_gdd.sh ${dir2}
echo "calc gsl"
. 013_gsl.sh ${dir2}
echo "calc summer days"
. 014_summerdays.sh ${dir2}
echo "calc temp means"
. 015_tempmeans.sh ${dir2}

