#!/bin/bash

#javac -cp "src/com/FileManager.class:/lustre/expphy/work/halla/sbs/bclary/extras/coatjava_4a.8.1/lib/clas/coat-libs-3.0-SNAPSHOT.jar" phiAnalysis.java

if [ "$1" = "PID"  ]
then
 	javac -cp /lustre/expphy/work/halla/sbs/bclary/extras/coatjava_4a.8.2/lib/clas/coat-libs-4.0-SNAPSHOT.jar -d . -sourcepath . phiAnalysis.java
	java -cp ".:src/com/*.class:/lustre/expphy/work/halla/sbs/bclary/extras/coatjava_4a.8.2/lib/clas/coat-libs-4.0-SNAPSHOT.jar" phiAnalysis $2
elif [ "$1" = "PHYS" ]
then
	javac -cp /lustre/expphy/work/halla/sbs/bclary/extras/coatjava_4a.8.2/lib/clas/coat-libs-4.0-SNAPSHOT.jar -d . -sourcepath . physicsAnalysis.java
	java -cp ".:src/com/*.class:/lustre/expphy/work/halla/sbs/bclary/extras/coatjava_4a.8.2/lib/clas/coat-libs-4.0-SNAPSHOT.jar" physicsAnalysis $2
fi
