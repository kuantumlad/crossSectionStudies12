#!/bin/bash

#scripts to submit jobs to batch farm


njobs=`ls -l /lustre/expphy/volatile/clas12/bclary/lund_files/elastic_electron/*.lund | wc -l`
file_prefix="clas12_2GeV_"

echo "Creating jobs for ${njobs} files"

i=0
rm *.jsub
while [ $i -lt $njobs ]
do

    outfile=sub_elast_${i}.jsub
    echo "JOBNAME: clas12ElasticClary${i}" >> $outfile
    echo "OS: centos7" >> $outfile
    echo "TRACK: simulation" >> $outfile
    echo "PROJECT: clas12" >> $outfile
    echo "MEMORY: 8096 MB" >> $outfile
    echo "INPUT_FILES: /lustre/expphy/volatile/clas12/bclary/lund_files/elastic_electron/clas12_2GeV_${i}.lund" >> $outfile
    echo "OTHER_FILES: /w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/elast_gen_cole/run-gemc.csh" >> $outfile
    echo "COMMAND: ./run-gemc.csh ${file_prefix}${i}.lund" >> $outfile
    echo "OUTPUT_DATA: out.hipo" >> $outfile
    echo "OUTPUT_TEMPLATE: /lustre/expphy/volatile/clas12/bclary/simulation_out/elastic/beam2/clas12_2GeV_${i}_tp1sm1.hipo" >> $outfile

    
i=$(($i+1))
echo $i
    
done

j=0
# submit jobs
while [ $j -lt $njobs ]
do
    echo "submitting job ${j} for file sub_elast_${j}.jsub "
    jsub sub_elast_${j}.jsub
    j=$((j+1))
done
