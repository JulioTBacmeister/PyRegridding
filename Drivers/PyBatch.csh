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


module load conda

conda activate npl-2022b

echo "Cruising .... "
./batch_drv.py
