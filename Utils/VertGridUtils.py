import numpy as np
import MyConstants as C
import xarray as xr

def Original_as_bs( px , pcross=0.086 , p58=None, gridalign=False ):
    #----------------------------------------------------------------------------------
    #  This returns hybrid A and B coefficients using what
    #  seems to be the accepted historical CESM method (ref?)
    #  for creating the vertical grid:
    #     1) Piecewise linear  A's and B's (in ref. press. or normalized ref. press.)
    #     2) Crossover to pure pressure coords somewhere around 100 hPa
    #  Inputs
    #     px:      normalized (1 to 0) ref. interface pressures, i.e., px = hyai+hybi 
    #               in NEW grid
    #     pcross:  Actual crossover press to pure p-coords.  This is here to account
    #               for a bit of a mess in CESM3 vertical grids at lowest layer
    #     p58:     normalized (1 to 0) ref. interface pressures for L58 CESM3 vertical
    #               grid. Needed if gridalign=True
    #     gridalign: Adjust V-grid to be consistent with CESM3 L58
    #----------------------------------------------------------------------------------
    
    m=( (1.0-0.0) / (1.0-pcross ) )
    Lx = len( px )

    if ( gridalign == True ):
        oo=np.where( p58 > pcross )
        pshift = pcross - p58[ oo[0][0]-1 ]
    else:
        pshift=0.0

    hybi_1 = m*(px -pcross +pshift)
    hyai_1 = px
    ##
    hybi_2 = np.where( hybi_1<0. , 0., hybi_1 )
    hyai_2 = px - hybi_2 # - m*pshift
    ###
    hyai_3 = np.where( hyai_2<0. , 0., hyai_2 )
    hybi_3 = px - hyai_3
    
    hyai,hybi = hyai_3, hybi_3
    return hyai,hybi