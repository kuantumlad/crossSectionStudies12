#!/bin/bash

# $1: dir in
# $2: run number
# $3: min file
# $4: max file

#  /work/clas12/clas12/data/trains/v2/

#Spring 2018:
#4013

#Fall 2018:
#5000, 5001, 5030, 5036, 5038, 5046, 5117

javac -cp /w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/coat/coatjava-5b.7.6/lib/clas/coat-libs-5.7.6-SNAPSHOT.jar  -d . -sourcepath . accpTest.java 
java -cp ".:/w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/coat/coatjava-5b.7.6/lib/clas/coat-libs-5.7.6-SNAPSHOT.jar" accpTest $1 $2 $3 $4
