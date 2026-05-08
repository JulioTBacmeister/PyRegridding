#!/usr/bin/env python
# Import packages 
import sys
import argparse as arg
# import modules in other directories
# sys.path.append('/glade/work/juliob/PyRegridding/Regridder/')
# sys.path.append('/glade/work/juliob/PyRegridding/Utils/')
sys.path.append('../Regridder/')
sys.path.append('../Utils/')


import importlib
import glob
import copy
#import time
import os 
import subprocess as sp

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri

try:
    import ESMF as E
except ImportError:
    import esmpy as E

import scripGen as SG
import esmfRegrid as erg
import scipy

# "ChatGPI version" --- 
import VertRegridFlexLL as vrg
print( "Using Flexible parallel/serial VertRegrid ")

import GridUtils as GrU
import MakePressures as MkP
import humiditycalcs as hum
import MyConstants as Con

# Reload local packages that are under
# development
importlib.reload( erg )
importlib.reload( vrg )
importlib.reload( SG )
importlib.reload( MkP )
importlib.reload( hum )
importlib.reload( GrU )
importlib.reload( Con )
#importlib.reload( Gv )



import glob
import subprocess
import os



def nccopy(BaseDir,ymdPat,Hgrid,Vgrid):

# Specify the directory path and the pattern to match files

#basedir = "/glade/derecho/scratch/juliob/ne240np4_IC_for_pg3"
basedir = "/glade/campaign/cgd/amp/juliob/QBO_nudging"

Vgrid ="L93"
Hgrid="ne30pg3"
rVgrid = f"/{Vgrid}/"
rVgrid_cdf5 = f"/{Vgrid}_cdf5/"
timetag =f"2002-*-*-*"

#directory_path = (
#    f"/glade/campaign/cgd/amp/juliob/ERA5/{Hgrid}/{Vgrid}/"  # '/path/to/directory'
#    )
directory_path = (
    f"{basedir}/{Hgrid}/{Vgrid}/"  # '/path/to/directory'
    )

#directory_path_cdf5 = directory_path.replace( "/L93/", "/L93_cdf5/" )
directory_path_cdf5 = directory_path.replace( rVgrid , rVgrid_cdf5 )

print( directory_path_cdf5 )

command = ["mkdir", "-p", directory_path_cdf5 ]
subprocess.run(command, check=True)

file_pattern = f"ERA5*{timetag}.nc"  # Example pattern to match NetCDF files, adjust as needed

# Use glob.glob() to generate a list of files matching the pattern in the directory
file_list = sorted(glob.glob(directory_path + "/" + file_pattern))


# print(file_list)

# Print the list of files
for file_i in file_list:
    #file_c = file_i.replace("/L93/", "/L93_cdf5/")
    file_c = file_i.replace( rVgrid , rVgrid_cdf5 )

    # Define the command
    command = ['ncdump', '-k', file_i]
    # Run the command and capture the output
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # Check if the command was successful
    if result.returncode == 0:
        # Print or use the output
        format_type = result.stdout.strip()
        print(f"NetCDF format of {file_i} is {format_type}")
    else:
    # Handle the error    
        print(f"Error: {result.stderr}")
        
    if (os.path.exists(file_c)==False):
        print(f"About to nccopy {file_i} .")
        # Construct the command to run
        command = ["nccopy", "-k", "cdf5", file_i, file_c]
        subprocess.run(command, check=True)
        print(f"nccopy {file_i} to {file_c} completed successfully.")
    else:
        print(f"{file_c} aready exists!!!")