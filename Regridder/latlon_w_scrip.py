#####
# Makes lat-lon vars from SCRIP file 
#####

import xarray as xr
import numpy as np
import pandas as pd


def latlon(scrip,gridHkey):
    
    S = xr.open_dataset( scrip )
    
    if (gridHkey == 'c'):
        lon = S.grid_center_lon.values
        lat = S.grid_center_lat.values
        
    if (gridHkey == 'yx'):
        lon = np.unique( S.grid_center_lon.values )
        lat = np.unique( S.grid_center_lat.values )
        
        
        
    return lat,lon