#!/bin/bash

#/volatile/clas12/rg-a/production/recon/calib/v14/005038/calib_clas_005038.evio.00240-00244.hipo 
current=240
current_plus=244


end=295

while [ $current -le $end ]
do
    # blocked into 5 files per file
    echo " converting files range ${current} to ${current_plus}"
    /work/clas12/tylern/software/hipo_tools/bin/hipo2root -t -tbt /volatile/clas12/rg-a/production/recon/calib/v14/005038/calib_clas_005038.evio.00${current}-00${current_plus}.hipo calib_clas_005038.evio.00${current}.root

    current=$((current+5))
    current_plus=$((current_plus+5))

done
