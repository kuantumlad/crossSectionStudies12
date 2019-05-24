#!/bin/csh


# need to get correct version of root that has h2root
source /site/12gev_phys/softenv.csh 2.1

set beamE = $1

# number of files to convert
set n = 0
while ( $n < 1000 )
    echo " converting file elast_gen_${beamE}GeV_${n}.rz to root "
    h2root elast_gen_${beamE}GeV_${n}.rz elast_gen_${beamE}GeV_${n}.root
    @ n++
end
