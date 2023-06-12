#/usr/bin/bash

python dataExtractorInterpolator.py --casedir "./Test1" --savedir "./Datasets" --datasetname "Test1.mat" --timestep 0.0002 --startstep 377579 --endstep 380525 --writeinterval 491 --Xmin -1.5 --Xmax 1.5 --Ymin -1.0 --Ymax 1.0 --positiondata False

python dataExtractorInterpolator.py --casedir "./Test2" --savedir "./Datasets" --datasetname "Test2.mat" --timestep 0.0002 --startstep 377579 --endstep 380525 --writeinterval 491 --Xmin -1.5 --Xmax 1.5 --Ymin -1.0 --Ymax 1.0 --positiondata False

python dataExtractorInterpolator.py --casedir "./Test3" --savedir "./Datasets" --datasetname "Test3.mat" --timestep 0.0002 --startstep 377579 --endstep 380525 --writeinterval 491 --Xmin -1.5 --Xmax 1.5 --Ymin -1.0 --Ymax 1.0 --positiondata False


