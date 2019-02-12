#!/bin/bash


target_dir=/volatile/clas/clase1/markov/HTCCCook/calib/
 

run=002391
num_files=`ls -l /volatile/clas/clase1/markov/HTCCCook/calib/*${run}*hipo | wc -l`
echo " Number of files " ${num_files}

f_min=$1
f_max=$2

min=$f_min

# naive error
if [ $f_max -gt $num_files ]
then
echo "Max files number exceeds number of files" 
fi

rm merge-files_${run}_${min}_${f_max}.sh

echo "#!/bin/bash" >> merge-files_${run}_${min}_${f_max}.sh 

printf "/work/clas12/bclary/CLAS12/validation/validationBrandon/coat/coatjava-5.7.4/bin/hipo-utils -merge -o /volatile/clas12/bclary/clas12RunData/run${run}/out_clas12_${run}.evio.${f_min}.${f_max}.hipo " >> merge-files_${run}_${min}_${f_max}.sh 

while [ $f_min -lt $f_max ]
do 

echo " File " $f_min 
f_min=$((f_min+1))

printf ${target_dir}"out_clas_"${run}".evio."${f_min}".hipo " >> merge-files_${run}_${min}_${f_max}.sh 

done 


chmod +x merge-files_${run}_${min}_${f_max}.sh 

./merge-files_${run}_${min}_${f_max}.sh 
