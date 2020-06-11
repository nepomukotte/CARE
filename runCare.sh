#!/bin/bash

if [[ "$LD_LIBRARY_PATH" =~ (^|:)"/usr/local/lib"(|/)(:|$) ]]; then
	:
else
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
fi

./CameraAndReadout -s 3242 -c SPB2/SPB2_CARE_Camera.cfg -of test -if /home/nepi/Share/POEMMA/SPB2/CARETRIGGERStudies/SimulatePointSourceInFP/SPB2PointSourceSims.root 
