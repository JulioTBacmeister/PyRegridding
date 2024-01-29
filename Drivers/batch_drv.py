#!/usr/bin/env python
# Import packages 
import os 
import subprocess as sp
import update_config as uc
import DrvRegrid as DR

#####################################################
# This function is called by PyBatch.csh, and
# then also submits PyBatch.csh after incrementing
# month in config file.
#####################################################

def main():

    theYear = 2003
    file_path = './config.yaml'  # Specify the path to your config file

    config = uc.read_config_yaml( file_path )
    print( config )
    
    # Add the regrid commands here:
    # ...
    #. ./DrvRegrid.py --year=2000 --month=$month --day=99 --hour=99 --Dst='ne30pg3' --DstVgrid='L93'
    
    DR.main( year=config['year'] , month=config['month'] , day=config['day'], hour=config['hour'] , Dst=config['Dst'] , DstVgrid=config['DstVgrid'] , Src='ERA5' )
    
    #------------------------------
    
    config = uc.increment_month( config )
    print( config )
    uc.write_config_yaml(file_path, config)
   
    if ((config['year']==theYear) and (config['month']<=12) and (config['Resubmit']>0) ):
        print(f" Resubmitting myself through PyBatch.csh  ")
        
        sp.run(f"qsub PyBatch.csh", 
               shell=True )
        print(f"PyBatch ... " )
        
    
    

if __name__ == "__main__":
    main()
