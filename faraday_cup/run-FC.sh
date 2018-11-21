#!/bin/bash

# $1: dir in
# $2: run number
# $3: min file
# $4: max file

javac -cp /w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/coat/coatjava-5b.7.1/lib/clas/coat-libs-5.7.1-SNAPSHOT.jar -d . -sourcepath . faradayCupCalculator.java 
java -cp ".:/w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/coat/coatjava-5b.7.1/lib/clas/coat-libs-5.7.1-SNAPSHOT.jar" faradayCupCalculator $1 $2 $3 $4