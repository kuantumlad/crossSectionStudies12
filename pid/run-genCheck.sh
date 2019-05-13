#!/bin/bash

javac -cp /w/hallb-scifs17exp/clas12/bclary/CLAS12/simulation/farm/coatjava-5b.7.8/lib/clas/coat-libs-5.7.8-SNAPSHOT.jar -d . -sourcepath . genCheck.java 
java -cp ".:/w/hallb-scifs17exp/clas12/bclary/CLAS12/simulation/farm/coatjava-5b.7.8/lib/clas/coat-libs-5.7.8-SNAPSHOT.jar" genCheck $1
