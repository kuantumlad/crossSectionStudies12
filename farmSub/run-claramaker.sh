#!/bin/bash

anaName=$1
num_clara_jobs=1

i=0

echo " >> Reconstruction for ${anaName} " 
rm toCook_*.clara

while [ $i -lt $num_clara_jobs ]
do
        
    echo "set inputDir /lustre/expphy/volatile/clas12/bclary/simulation_out/elastic/beam2/" >> toCook_${i}.clara
    echo "set outputDir /work/clas12/bclary/CLAS12/electron_studies/elastic/sim/recon/"  >> toCook_${i}.clara
    echo "set description clary" >> toCook_${i}.clara 
    echo "set farm.track reconstruction" >> toCook_${i}.clara 
    echo "set farm.scaling 15" >> toCook_${i}.clara 
    echo "set farm.cpu 70" >> toCook_${i}.clara 
    echo "set farm.memory 60" >> toCook_${i}.clara 
    echo "set threads 64"  >> toCook_${i}.clara 
    echo "set fileList file_list.list"  >> toCook_${i}.clara 
    echo "run farm"  >> toCook_${i}.clara 
  
    /work/clas12/bclary/CLAS12/validation/validationBrandon/clara-5c.7.0/bin/clara-shell toCook_${i}.clara 

    i=$((i+1))

done
