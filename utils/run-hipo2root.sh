#!/bin/bash

inDir=$1


#tot_files=`ls /work/clas12/bclary/CLAS12/electron_studies/elastic/sim/recon/field_tp1sm1/*.hipo | wc -l`
#path_to_files="/work/clas12/bclary/CLAS12/electron_studies/elastic/sim/recon/field_tp1sm1/"

#inclusive path  n file name
#/volatile/clas12/bclary/simulation_out/inclusive/7GeV/norad/out_lund_files_n462.dat.hipo

tot_files=1000 #`ls ${inDir}*.hipo | wc -l`
path_to_files=${inDir}


echo " Converting " ${tot_files} " hipo files to root files " 

i=1
while [ $i -lt ${tot_files} ]
do

    echo " convert file ${path_to_files}${i}.dat.hipo ---> ${path_to_files}${i}.root  "
    #/work/clas12/bclary/CLAS12/validation/validationBrandon/hipo2root/hipo2root ${path_to_files}out_clas12_2GeV_${i}_tp1sm1.hipo ${path_to_files}out_clas12_2GeV_${i}_tp1sm1.root 
    /work/clas12/bclary/CLAS12/validation/validationBrandon/hipo2root/hipo2root ${path_to_files}${i}.lund.hipo ${path_to_files}${i}.root 
    i=$((i+1))

done

echo " Done converting files"
