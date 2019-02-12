#!/bin/csh -f
source /group/clas12/gemc/environment.csh 4a.2.2
setenv line $1
echo $line
~/hipo2root/hipo2root /work/clas12/clas12/data/calib/cooked_5bp6p1/$line $line.root
/u/home/markov/filtering/hello2GeV $line.root filtered$line.root
mv $line.root /lustre/expphy/volatile/clas/clase1/markov/12GeV/inclusive/inclusiveDiffPol/farm/filesClara/10.6GeV/InclusiveElastic/data/6GeV
mv filtered$line.root /lustre/expphy/volatile/clas/clase1/markov/12GeV/inclusive/inclusiveDiffPol/farm/filesClara/10.6GeV/InclusiveElastic/data/6GeV
