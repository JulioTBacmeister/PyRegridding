#!/bin/csh
#PBS -N sefvpy2

### Charging account
#PBS -A P93300642 
### Request one chunk of resources with N CPU and M GB of memory
#PBS -l select=1:ncpus=1:mem=64GB
### 
#PBS -l walltime=12:00:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe



# This job's working directory
#


module load conda

conda activate npl-2023b


# Set up to regrid one of the lates eighth degree runs to a regular lat-lon 1/8x1/8 grid.
# All 3 of the high frequency history streams are regridded h1i, h2i, and h3i
#----------------------------------------------------------------------------------------
./SE_x_FV.py --case=c153_topfix_ne240pg3_FMTHIST_xic_x02 --hsPat=cam.h1i --Src=ne240pg3 --Dst=latlonOxO --ymdPat='2004-07-*'
./SE_x_FV.py --case=c153_topfix_ne240pg3_FMTHIST_xic_x02 --hsPat=cam.h2i --Src=ne240pg3 --Dst=latlonOxO --ymdPat='2004-07-*'
./SE_x_FV.py --case=c153_topfix_ne240pg3_FMTHIST_xic_x02 --hsPat=cam.h3i --Src=ne240pg3 --Dst=latlonOxO --ymdPat='2004-07-*'



