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

#./SE_x_FV.py --case=POLARRES_FMT_x01_E --hsPat=cam.h1 --Src=POLARRES
#./SE_x_FV.py --case=cDEV_ne120pg3_FMTHIST_DYAMOND_x01_e02 --hsPat=cam.h1i --Src=ne120pg3 --Dst=fvQxQ
#./SE_x_FV.py --case=cDEV_ne240pg3_FMTHIST_aicn_x02 --hsPat=cam.h0a --Src=ne240pg3 --Dst=latlonOxO
#./SE_x_FV.py --case=cDEV_ne240pg3_FMTHIST_aicn_x02 --hsPat=cam.h1i --Src=ne240pg3 --Dst=fv1x1 --DstSubDir=regridded_1x1 --AllConservative
#./SE_x_FV.py --case=c153_ne30pg3_FMTHIST_x02 --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1

#./SE_x_FV.py --case=c153_ne240pg3_FMTHIST_aicn_x01 --hsPat=cam.h1i --Src=ne240pg3 --Dst=latlonOxO --ymdPat='2004-09-*'

#./SE_x_FV.py --case=c153_topfix_ne120pg3_FMTHIST_xic_x02 --hsPat=cam.h0a --Src=ne120pg3 --Dst=fv1x1 --ymdPat='*'
#./SE_x_FV.py --case=c153_topfix_ne240pg3_FMTHIST_xic_x02 --hsPat=cam.h0a --Src=ne240pg3 --Dst=fv1x1 --ymdPat='*'

#./SE_x_FV.py --case=c153_topfix_ne240pg3_FMTHIST_xic_x02 --hsPat=cam.h1i --Src=ne240pg3 --Dst=latlonOxO --ymdPat='2004-07-*'
#./SE_x_FV.py --case=c153_topfix_ne240pg3_FMTHIST_xic_x02 --hsPat=cam.h2i --Src=ne240pg3 --Dst=latlonOxO --ymdPat='2004-07-*'
#./SE_x_FV.py --case=c153_topfix_ne240pg3_FMTHIST_xic_x02 --hsPat=cam.h3i --Src=ne240pg3 --Dst=latlonOxO --ymdPat='2004-07-*'

#./SE_x_FV.py --case=c153_topfix_ne240pg3_FMTHIST_L56_86km_aicn_x01 --hsPat=cam.h1i --Src=ne240pg3 --Dst=latlonOxO --ymdPat='*'


#./SE_x_FV.py --case=c153_ne240pg3_FMTHIST_aicn_x01 --hsPat=cam.h0a --Src=ne240pg3 --Dst=latlonOxO --ymdPat='*'

#./SE_x_FV.py --case=c153_ne30pg3_FMTHIST_x02 --hsPat=cam.h1i --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'

#./SE_x_FV.py --case=c160mm_ne30pg3_FMTndg_cciasi_x05 --hsPat=cam.h1i.TAU_uniform --Src=ne30pg3 --Dst=CCIASI  --DstSubDir=CCIASI --ymdPat='*'
#./SE_x_FV.py --case=c160mm_ne30pg3_FMTndg_cciasi_x05 --hsPat=cam.h1i.TAU_vort500 --Src=ne30pg3 --Dst=CCIASI  --DstSubDir=CCIASI --ymdPat='*'

#./SE_x_FV.py --case=c64_005mm_ne30pg3_FMTHIST_CTL  --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'

#./SE_x_FV.py --case=c64_005_ne30pg3_FMTHIST_rdgres_x01 --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'
#./SE_x_FV.py --case=c64_005_ne30pg3_FMTHIST_CTL_x01 --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'
#./SE_x_FV.py --case=c64_005_ne30pg3_FMTHIST_C5phys_x06 --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'

#./SE_x_FV.py --case=c64_gwR2_ne30pg3_FMTHIST_topfix_x01   --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'
./SE_x_FV.py --case=c64_gwR2_ne30pg3_FMTHIST_topfix_oldCTL_qbotune_x03   --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'

#./SE_x_FV.py --case=cDEV_ne30pg3_FMTHIST_C5_x01 --hsPat=cam.h0a --Src=ne30pg3 --Dst=fv1x1 --ymdPat='*'



