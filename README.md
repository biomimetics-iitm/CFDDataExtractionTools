# dataExtractionInterpolation-CFD
This repository shall house the data extraction, inteprolation and other miscellaneous pre processing scripts used for (a) Data storage (b) data post processing and (c) data generation in ML work loads. 


## About the code

dataExtractorInterpolator.py

Requirements:
 - python version > 3.x
 - numpy
 - scipy
 
What does it do:
 - Loads AllInOne.dat files of IBM solver - "dt" needs to be specified. 
 - Select and extract a region of interest through user specification
   - Xmin, Xmax, Ymin, Ymax being the bounds of the bounding box
   - Specify the start and end time step with write interval as used in the case.
 - Save the whole data matrix into a mat file. 
 - Can be saved into h5 files as well (Requires h5py python package)
