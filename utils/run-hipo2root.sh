#!/bin/bash

inDir=$1

tot_files=66 #`ls ${inDir}*.hipo | wc -l`
path_to_files=${inDir}


echo " Converting " ${tot_files} " hipo files to root files " 

i=1
while [ $i -lt ${tot_files} ]
do

    echo " convert file ${path_to_files}${i}.dat.hipo ---> ${path_to_files}${i}.root  "
    /work/clas12/tylern/software/hipo_tools/bin/hipo2root -t -c ${path_to_files}${i}.hipo ${path_to_files}${i}.root 
    i=$((i+1))

done

echo " Done converting files"
