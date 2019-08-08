#!/bin/bash

javac -cp /w/hallb-scifs17exp/clas12/bclary/CLAS12/simulation/farm/coatjava-6.3.1/lib/clas/coat-libs-6.3.1-SNAPSHOT.jar -d . -sourcepath . genCheck.java 
java -cp ".:/w/hallb-scifs17exp/clas12/bclary/CLAS12/simulation/farm/coatjava-6.3.1/lib/clas/coat-libs-6.3.1-SNAPSHOT.jar" genCheck $1
