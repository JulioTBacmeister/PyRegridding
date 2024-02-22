#!/bin/csh
#PBS -N sefvpy

### Charging account
#PBS -A P93300642 
### Request one chunk of resources with N CPU and M GB of memory
#PBS -l select=1:ncpus=16:mem=256GB
### 
#PBS -l walltime=8:00:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe



# This job's working directory
#


module load conda

conda activate npl-2023b

./SE_x_FV.py --case=fmthist_MM_control_TEM --hsPat=cam.h1
./SE_x_FV.py --case=fmthist_MM_control_TEM --hsPat=cam.h0
