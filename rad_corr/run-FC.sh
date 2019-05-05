#!/bin/bash

# $1: dir in
# $2: run number
# $3: min file
# $4: max file

javac -cp /w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/coat/coatjava-5b.7.6/lib/clas/coat-libs-5.7.6-SNAPSHOT.jar  -d . -sourcepath . radCorr.java 
java -cp ".:/w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/coat/coatjava-5b.7.6/lib/clas/coat-libs-5.7.6-SNAPSHOT.jar" radCorr $1
