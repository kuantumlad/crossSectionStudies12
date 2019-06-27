#!/bin/bash

start=$1
end=$2
inputDir=$3

source /site/12gev_phys/1.2/Linux_CentOS6.5-x86_64-gcc4.4.7/root/5.34.21/root/bin/thisroot.csh

# elast_gen_2GeV_


while [ $start -lt $end ]
do

    echo " processing file ${start} "   
    /site/12gev_phys/1.2/Linux_CentOS6.5-x86_64-gcc4.4.7/root/5.34.21/root/bin/h2root ${inputDir}_${start}.rz ${inputDir}_${start}.root
    start=$[${start}+1]

done
