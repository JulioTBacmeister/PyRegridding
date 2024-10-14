#!/usr/bin/env python
# Import packages 
import os

from mpi4py import MPI


import subprocess as sp
import update_config as uc
import DrvRegrid as DR

######################################################################
# This function is called by PyBatch_ERA5regrid.csh, and
# then may also resubmit PyBatch_ERA5regrid.csh after incrementing
# month and decrementing Resubmit in config_ERA5regrid file.
#####################################################################

def main():

    # This what the yaml file should look like:
    """
    Dst: ne240np4
    DstVgrid: L93
    IC_for_pg: true
    Resubmit: 0
    StepBy: hour
    StepN: 6
    TheProcYear: 2004
    day: 1
    hour: 0
    month: 7
    year: 2004
    """
    
    file_path = './config_ERA5regrid.yaml'  # Specify the path to your config file

    config = uc.read_config_yaml( file_path )
    print( config )
    theYear = config['TheProcYear']
    # Destination grid
    Dst=config['Dst']
    
    if (Dst == 'ne480np4'):
        BestRegridMethod = 'CONSERVE_2ND' 
    elif (Dst == 'ne480pg3'):
        BestRegridMethod = 'CONSERVE_2ND' 
    elif (Dst == 'ne240np4'):
        BestRegridMethod = 'CONSERVE_2ND' 
    elif (Dst == 'ne240pg3'):
        BestRegridMethod = 'CONSERVE_2ND' 
    elif (Dst == 'ne120np4'):
        BestRegridMethod = 'CONSERVE_2ND' 
    elif (Dst == 'ne120pg3'):
        BestRegridMethod = 'CONSERVE_2ND' 
    elif (Dst == 'ne30np4'):
        BestRegridMethod = 'CONSERVE_2ND' 
    elif (Dst == 'ne30pg3'):
        BestRegridMethod = 'CONSERVE_2ND' 
    else:
        BestRegridMethod = 'CONSERVE' 

    print(f' "I" have decided that the best regrid method is {BestRegridMethod}')
    
    # Add the regrid commands here:
    # ...
    #. ./DrvRegrid.py --year=2000 --month=$month --day=99 --hour=99 --Dst='ne30pg3' --DstVgrid='L93'
    
    DR.main( year=config['year'] , month=config['month'] , day=config['day'], hour=config['hour'] , Dst=config['Dst'] , DstVgrid=config['DstVgrid'] , IC_for_pg=config['IC_for_pg'], Src='ERA5' , RegridMethod=BestRegridMethod )
    
    #------------------------------
    if (config['StepBy'].lower() == 'day'):
        config = uc.increment_day( config ) #, NoLeapYear=True )
    if (config['StepBy'].lower() == 'month'):
        config = uc.increment_month( config ) #, NoLeapYear=True )
    if (config['StepBy'].lower() == 'hour'):
        config = uc.increment_hours( config , nhours=config['StepN']) #, NoLeapYear=True )

    config = uc.decrement_Resubmit( config )
    print( config )
    uc.write_config_yaml(file_path, config)
   
    if (( (config['year']==theYear) or (theYear<0) ) and (config['month']<=12) and (config['Resubmit']>=0) ):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        if (rank == 0):
            print(f" Resubmitting myself through PyBatch.csh  ")
            
            sp.run(f"qsub PyBatch_ERA5regrid.csh", 
                   shell=True )
            print(f"PyBatch ... " )
        
    
    

if __name__ == "__main__":
    main()
