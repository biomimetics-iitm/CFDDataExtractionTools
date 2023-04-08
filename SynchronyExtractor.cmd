#!/bin/bash
#PBS -N "Extractor"
#PBS -e errorfile.err
#PBS -o logfile.log
#PBS -q workq
#PBS -l select=1:host=node1:ncpus=8
#PBS -V
cd $PBS_O_WORKDIR

module load conda
source activate tf-gpu

/home/sunetra/.conda/envs/tf-gpu/bin/python3 <path to the parent folder folder>dataExtractorInterpolator.py > outputlog.txt
