#!/bin/csh -f
source /group/clas12/gemc/environment.csh 4a.2.2
setenv line $1
echo $line
~/hipo2root/hipo2root /lustre/expphy/volatile/clas/clase1/markov/2GeV/inclusive/elastgen/farm/files3Particles/diffFields/2GeV/5.6.1/$line $line.root
./hello2GeV $line.root filtered$line.root
mv $line.root /lustre/expphy/volatile/clas/clase1/markov/2GeV/inclusive/elastgen/farm/files3Particles/diffFields/2GeV/5.6.1/
mv filtered$line.root /lustre/expphy/volatile/clas/clase1/markov/2GeV/inclusive/elastgen/farm/files3Particles/diffFields/2GeV/5.6.1/
hadd /lustre/expphy/volatile/clas/clase1/markov/2GeV/inclusive/elastgen/farm/files3Particles/diffFields/2GeV/5.6.1/filtered*.root /lustre/expphy/volatile/clas/clase1/markov/2GeV/inclusive/elastgen/farm/files3Particles/diffFields/2GeV/5.6.1/filteredT-1S1.root
