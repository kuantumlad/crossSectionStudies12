#!/bin/csh -f
setenv run $1
echo $run
foreach line("`cat run.dat`")
setenv outFile jsub/out$line.jsub
echo PROJECT: clas12 >  $outFile
echo TRACK: simulation >>   $outFile
echo JOBNAME: $line>>  $outFile
echo MEMORY: 8000 MB>>  $outFile
echo CPU: 1 >> $outFile
echo OS: centos7>>  $outFile
echo MAIL: markov@jlab.org>>  $outFile
echo SINGLE_JOB : TRUE>>  $outFile
echo OTHER_FILES: /u/home/markov/filtering/convertSimple.sh>>  $outFile
echo COMMAND: ./convertSimple.sh $line >>$outFile
jsub $outFile
end

