#!/bin/bash


# /work/clas12/clas12/data/calib/cooked_5bp7p1/out_clas_002391.evio.1.hipo
dir_path_in=/work/clas12/clas12/data/calib/cooked_5bp7p1
tot_files=200 #`ls ${dir_path_in}/out_clas_002391.*.hipo | wc -l`
#tot_files=`ls /work/clas12/bclary/CLAS12/electron_studies/elastic/sim/recon/out_clas12*.hipo | wc -l`
dir_path=/work/clas12/bclary/CLAS12/electron_studies/elastic/data
converter_path=/w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/hipo2root

echo " Total number of to convert $tot_files "
file_number=0

while [ ${file_number} -lt ${tot_files} ]
do

    echo " converting file $file_number "
    ${converter_path}/hipo2root ${dir_path_in}/out_clas_002391.evio.${file_number}.hipo ${dir_path}/out_clas12_002391.${file_number}.root
    file_number=$((file_number+1))

done
