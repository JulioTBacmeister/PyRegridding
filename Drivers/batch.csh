#!/bin/csh
#PBS -N ERA5proc

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
echo `date`
if ( $?PBS_O_WORKDIR == 1 ) then
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
endif



module load conda

conda activate npl-2022b

source config.txt

#./DrvRegrid.py --year=2021 --month=6 --day=99 --hour=99 --Dst='Arctic' --DstVgrid='L32'
./DrvRegrid.py --year=2000 --month=$month --day=99 --hour=99 --Dst='ne30pg3' --DstVgrid='L93'
