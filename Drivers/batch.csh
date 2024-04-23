#!/bin/csh
#PBS -N sefvpy2

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

#./SE_x_FV.py --case=POLARRES_FMT_x01_E --hsPat=cam.h1 --Src=POLARRES
#./SE_x_FV.py --case=cDEV_ne120pg3_FMTHIST_DYAMOND_x01_e02 --hsPat=cam.h1i --Src=ne120pg3 --Dst=fvQxQ
#./SE_x_FV.py --case=cDEV_ne240pg3_FMTHIST_aicn_x02 --hsPat=cam.h0a --Src=ne240pg3 --Dst=latlonOxO

./SE_x_FV.py --case=cDEV_ne240pg3_FMTHIST_aicn_x02 --hsPat=cam.h1i --Src=ne240pg3 --Dst=fv1x1 --DstSubDir=regridded_1x1 --AllConservative


