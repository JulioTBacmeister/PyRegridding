#!/usr/bin/env python
# Import packages 
import os 
import subprocess as sp
import update_config as uc
import DrvRegrid as DR

######################################################################
# This function is called by PyBatch_ERA5regrid.csh, and
# then may also resubmit PyBatch_ERA5regrid.csh after incrementing
# month and decrementing Resubmit in config_ERA5regrid file.
#####################################################################

def main():

    file_path = './config_ERA5regrid.yaml'  # Specify the path to your config file

    config = uc.read_config_yaml( file_path )
    print( config )
    theYear = config['TheProcYear']
    
    # Add the regrid commands here:
    # ...
    #. ./DrvRegrid.py --year=2000 --month=$month --day=99 --hour=99 --Dst='ne30pg3' --DstVgrid='L93'
    
    DR.main( year=config['year'] , month=config['month'] , day=config['day'], hour=config['hour'] , Dst=config['Dst'] , DstVgrid=config['DstVgrid'] , Src='ERA5' )
    
    #------------------------------
    if (config['StepBy'].lower() == 'day'):
        config = uc.increment_day( config ) #, NoLeapYear=True )
    if (config['StepBy'].lower() == 'month'):
        config = uc.increment_month( config ) #, NoLeapYear=True )
    config = uc.decrement_Resubmit( config )
    print( config )
    uc.write_config_yaml(file_path, config)
   
    if ((config['year']==theYear) and (config['month']<=12) and (config['Resubmit']>0) ):
        print(f" Resubmitting myself through PyBatch.csh  ")
        
        sp.run(f"qsub PyBatch_ERA5regrid.csh", 
               shell=True )
        print(f"PyBatch ... " )
        
    
    

if __name__ == "__main__":
    main()
