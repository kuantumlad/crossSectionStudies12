#!/bin/tcsh 

rm elast_gen.* 
cp ../elast_gen . 
#elast_gen < ../inputFiles/elas_2036.inp 
elast_gen < ../inputFiles/e1f.inp 
