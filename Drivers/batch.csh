#!/bin/csh
#PBS -N sefvpy2

### Charging account
#PBS -A P93300042 
### Request one chunk of resources with N CPU and M GB of memory
#PBS -l select=1:ncpus=1:mem=128GB
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


#----------------------------------------------------------------------------------------
./SE_x_FV.py --case=f.e13.FAMIPC5.ne120_ne120_mt12.DECK_1950TR_RCP85_CDR.006 --hsPat=cam.h0 --Src=ne120np4 --Dst=latlonQxQ --ymdPat='*'
./SE_x_FV.py --case=f.e13.FAMIPC5.ne120_ne120_mt12.DECK_1950TR_CDR_hourly.006 --hsPat=cam.h0 --Src=ne120np4 --Dst=latlonQxQ --ymdPat='*'

#./SE_x_FV.py --case=CTL_fmt_c6477 --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'
#./SE_x_FV.py --case=CTL_fmt_c6477 --hsPat=cam.h1a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'

#./SE_x_FV.py --case=tausrf150_fmt_c6477 --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'
#./SE_x_FV.py --case=tausrf150_fmt_c6477 --hsPat=cam.h1a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'

#./SE_x_FV.py --case=tausrf200_fmt_c6477 --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'
#./SE_x_FV.py --case=tausrf200_fmt_c6477 --hsPat=cam.h1a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'

