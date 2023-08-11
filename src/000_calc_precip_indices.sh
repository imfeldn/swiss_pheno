#!/bin/bash

set -o errexit
set -o pipefail

dir2="../data/climate_indices/"

echo "calc cwd,cdd"
. 001_cwd_cdd.sh ${dir2}
echo "calc prcptot"
. 002_prcptot.sh ${dir2}
echo "calc rx5day"
. 003_rx5day.sh ${dir2}
echo "calc wd"
. 004_wd.sh ${dir2}
echo "calc snowdays"
. 005_snowdays.sh ${dir2}

