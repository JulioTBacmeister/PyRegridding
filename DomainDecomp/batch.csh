#!/bin/csh
#PBS -N pllel

### Charging account
#PBS -A P93300642 
### Request one chunk of resources with N CPU and M GB of memory
#PBS -l select=1:ncpus=4:mem=256GB
### 
#PBS -l walltime=0:30:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe


mpirun -np 4 ./stoopyExample.py
