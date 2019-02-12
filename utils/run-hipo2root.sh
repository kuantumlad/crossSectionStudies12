#!/bin/bash


tot_files=`ls /work/clas12/bclary/CLAS12/electron_studies/elastic/sim/recon/field_tp1sm1/*.hipo | wc -l`
path_to_files="/work/clas12/bclary/CLAS12/electron_studies/elastic/sim/recon/field_tp1sm1/"


echo " Converting " ${tot_files} " hipo files to root files " 

i=1
while [ $i -lt ${tot_files} ]
do

    echo " convert file " $i
    /work/clas12/bclary/CLAS12/validation/validationBrandon/hipo2root/hipo2root ${path_to_files}out_clas12_2GeV_${i}_tp1sm1.hipo ${path_to_files}out_clas12_2GeV_${i}_tp1sm1.root 
    i=$((i+1))

done

echo " Done converting files"
