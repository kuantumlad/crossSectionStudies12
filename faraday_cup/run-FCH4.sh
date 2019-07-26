#!/bin/bash

# $1: dir in
# $2: run number
# $3: min file
# $4: max file

#/volatile/clas12/rg-a/production/recon/calib/v14/005038/calib_clas_005038.evio.01050-01054.hipo

javac -cp /work/clas12/bclary/extras/coatjava-6b.3.0/lib/clas/coat-libs-6.3.0-SNAPSHOT.jar  -d . -sourcepath . faradayCupCalculatorH4.java 
java -cp ".:/work/clas12/bclary/extras/coatjava-6b.3.0/lib/clas/coat-libs-6.3.0-SNAPSHOT.jar" faradayCupCalculatorH4 $1 $2 $3 $4
