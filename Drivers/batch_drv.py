#!/usr/bin/env python
# Import packages 
import os 
import subprocess as sp
import update_config as uc
import DrvRegrid as DR

def main():

    theYear = 2002
    file_path = './config.txt'  # Specify the path to your config file

    config = uc.read_config( file_path )
    print( config )
    
    # Add the regrid commands here:
    # ...
    #. ./DrvRegrid.py --year=2000 --month=$month --day=99 --hour=99 --Dst='ne30pg3' --DstVgrid='L93'
    
    DR.main( year=config['year'] , month=config['month'] , day=99, hour=99 , Dst='ne30pg3' , DstVgrid='L93' , Src='ERA5' )
    
    #------------------------------
    
    config = uc.increment_month( config )
    print( config )
    uc.write_config(file_path, config)
   
    if ((config['year']==theYear) and (config['month']<=12)):
        print("YAAA AAAAaaaa PP ")
        
        sp.run(f"qsub PyBatch.csh", 
               shell=True )
        print(f"PyBatch ... " )
        
    
    

if __name__ == "__main__":
    main()
