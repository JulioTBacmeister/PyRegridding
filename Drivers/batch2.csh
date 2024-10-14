#!/bin/csh
#PBS -N sefvpy4

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


./SE_x_FV.py --case=c153_ne30pg3_FMTHIST_x05 --hsPat=cam.h1i --Src=ne30pg3 --Dst=fv1x1 --DstSubDir=regrid1x1 --ymdPat='2004-06-*'
./SE_x_FV.py --case=c153_topfix_ne120pg3_FMTHIST_xic_x02 --hsPat=cam.h1i --Src=ne120pg3 --Dst=fv1x1 --DstSubDir=regrid1x1 --ymdPat='2004-06-*'
./SE_x_FV.py --case=c153_topfix_ne240pg3_FMTHIST_xic_x02 --hsPat=cam.h1i --Src=ne240pg3 --Dst=fv1x1 --DstSubDir=regrid1x1 --ymdPat='2004-06-*'


