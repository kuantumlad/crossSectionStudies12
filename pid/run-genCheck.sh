#!/bin/bash

javac -cp /w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/coat/coatjava-5c.7.0/lib/clas/coat-libs-5.7.0-SNAPSHOT.jar -d . -sourcepath . genCheck.java 
java -cp ".:/w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/coat/coatjava-5c.7.0/lib/clas/coat-libs-5.7.0-SNAPSHOT.jar" genCheck $1
