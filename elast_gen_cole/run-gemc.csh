#!/bin/csh

# read in inputs
#set i=$1
set input_file=$1

setenv JLAB_ROOT /site/12gev_phys
source $JLAB_ROOT/softenv.csh 2.3

echo " Copying gcard over"
cp /w/hallb-scifs17exp/clas12/bclary/CLAS12/electron_studies/elast_gen_cole/clas12.gcard .
cp -r /w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/coat/coatjava-5c.7.0 .

echo " Getting events to process "
set NEV=`awk 'END {print NR}' ${input_file}`

echo "NUMBER OF EVENTS TO PROCESS " ${NEV}
echo "File to process ${input_file}"

# outbending t +1 and s -1
# inbending t - 1 and s -1
/group/clas12/gemc/4a.2.5/source/gemc clas12.gcard -N=${NEV} -USE_GUI=0 -RUNNO=11 -INPUT_GEN_FILE="LUND, ${input_file}" -SCALE_FIELD="TorusSymmetric, 1" -SCALE_FIELD="clas12-newSolenoid, -1"

echo " Decoding file now "

coatjava-5c.7.0/bin/evio2hipo  -r 11 -t 1 -s -1 out.ev -o out.hipo

echo "done"