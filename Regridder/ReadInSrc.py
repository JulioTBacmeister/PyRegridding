#!/usr/bin/env python
# Import packages 
import sys

# This contains all the data
#----------------------------
import GlobalVarClass
from GlobalVarClass import Gv


import xarray as xr
import numpy as np
import pandas as pd

try:
    import ESMF as E
except ImportError:
    import esmpy as E

import importlib
import copy
import time

import esmfRegrid as erg

#import cfgrib


import dask
import dask.array as da


# import modules in other directories
sys.path.append('../Utils/')
import GridUtils as GrU
import MakePressures as MkP
import humiditycalcs as hum
import MyConstants as Con

# Reload local packages that are under
# development
importlib.reload( erg )
importlib.reload( MkP )
importlib.reload( hum )
importlib.reload( GrU )
importlib.reload( Con )

# Physical Constants
Rdry = Con.Rdry() # 

#-------------------------------------------------------------
#  Naming conventions
#-------------------------------------------------------------
# aaa_{CAM,ERA} 
# Indicates the immediate provenance of a variable, e.g.,
#      phis_CAM ==> phis from CAM on CAM's grid
#      phis_ERA ==> phis from ERA on ERA's grid
# lower case 'phis' indicates this is ndarray-like 
#
# aaa_{CAM,ERA}_x{ERA,CAM}
# Indicates variable has been remapped horizontall. So, e.g.
#      phis_ERA_xCAM ==> ERA phis remapped horizontally to the CAM grid 
#
# aaa_{CAM,ERA}_xz{ERA,CAM}
# Indicates variable has been remapped horizontally AND vertically. So, e.g.
#      te_ERA_xzCAM ==> ERA temperature remapped horizontally to the CAM horizontal grid 
#                       and then also vertically interpoated to the CAM vertical grid
#-------------------------------------------------------------

# Define a function that loads a NetCDF file and returns an xarray dataset
@dask.delayed
def load_file(path):
    ds = xr.open_mfdataset(path ,  data_vars='different', coords='different' )
    return ds


def get_ERA5( year=2022, month=11, day=1, hour0=99):

    try:
        Gv.MySrc
        if ( Gv.MySrc != 'ERA5'):
            print( "You shouldnt be here - ABORT")
            rcode=-1
            return rcode
    except NameError:
        print( 'go on ' )

    tic_overall = time.perf_counter()

    monStr=str( year ).zfill(4)+str(month).zfill(2)
    # CAM history style yyyy-mm-dd string
    ymdStr=str( year ).zfill(4) + '-' + str(month).zfill(2) + '-' + str(day).zfill(2)

    if ( hour0 != 99 ):
        hour1=hour0+5
        ymdh0=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+str(hour0).zfill(2)
        ymdh1=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+str(hour1).zfill(2)
        ymdh=ymdh0+'_'+ymdh1
    else:
        ymdh0=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+"*"
        ymdh=ymdh0
        
    print( "Time tags for ERA5 files ...")
    print(monStr)
    print(ymdh) 

    era5dir = "/glade/campaign/collections/rda/data/ds633.6/e5.oper.an.ml/"
    wrkdir=era5dir+monStr+"/"

    print(Gv.dstTZHkey)

    #Define all file names for later use in dask function
    #-----------------------------------------------------
    spfile= wrkdir + 'e5.oper.an.ml.128_134_sp.regn320sc.'+ymdh+'.nc'
    tfile = wrkdir + 'e5.oper.an.ml.0_5_0_0_0_t.regn320sc.'+ymdh+'.nc'
    qfile = wrkdir + 'e5.oper.an.ml.0_5_0_1_0_q.regn320sc.'+ymdh+'.nc'
    ufile = wrkdir + 'e5.oper.an.ml.0_5_0_2_2_u.regn320uv.'+ymdh+'.nc'
    vfile = wrkdir + 'e5.oper.an.ml.0_5_0_2_3_v.regn320uv.'+ymdh+'.nc'
    wfile = wrkdir + 'e5.oper.an.ml.0_5_0_2_8_w.regn320sc.'+ymdh+'.nc'
    all_ERA_files = [ spfile , tfile, qfile, ufile, vfile, wfile ]

    """
    print( "Using DASK " )
    # Create a list of delayed objects, one for each file
    delayed_datasets = [load_file(path) for path in  all_ERA_files  ]
    # Use dask.compute to load all of the datasets in parallel
    datasets  = dask.compute(*delayed_datasets)
    dsPS_ERA  = datasets[0] 
    dsT_ERA   = datasets[1] 
    dsQ_ERA   = datasets[2] 
    dsU_ERA   = datasets[3] 
    dsV_ERA   = datasets[4] 
    dsW_ERA   = datasets[5] 

    """
    #Serial read
    #--------------
    print( 'Serial read of data ' )
    tic = time.perf_counter()
    dsPS_ERA  = xr.open_mfdataset( spfile, data_vars='different', coords='different' )
    dsT_ERA   = xr.open_mfdataset( tfile , data_vars='different', coords='different')
    dsQ_ERA   = xr.open_mfdataset( qfile,  data_vars='different', coords='different' )
    dsU_ERA   = xr.open_mfdataset( ufile , data_vars='different', coords='different')
    dsV_ERA   = xr.open_mfdataset( vfile , data_vars='different', coords='different')
    dsW_ERA   = xr.open_mfdataset( wfile , data_vars='different', coords='different')
    toc = time.perf_counter()
    print( all_ERA_files )
    pTime = f"Reading data vars for {ymdh} took  {toc - tic_overall:0.4f} seconds"
    print(pTime)

   

    tic = time.perf_counter()
    Gv.ps_ERA = dsPS_ERA['SP'].values
    Gv.te_ERA = dsT_ERA['T'].values
    Gv.q_ERA  = dsQ_ERA['Q'].values
    Gv.u_ERA  = dsU_ERA['U'].values
    Gv.v_ERA  = dsV_ERA['V'].values
    Gv.w_ERA  = dsW_ERA['W'].values

    Gv.lon_ERA = dsT_ERA['longitude'].values
    Gv.lat_ERA = dsT_ERA['latitude'].values
    toc = time.perf_counter()
    pTime = f"Extracting values from Xarrays took {toc - tic_overall:0.4f} seconds"
    print(pTime)
   
    #---------------------------
    # Get shape of ERA data
    #-----------------------------
    nt,nL,nx,ny = np.shape( Gv.te_ERA )

    #    Have a look at ERA global mean surface pressures
    for n in np.arange(nt):
        globPS = np.sum( Gv.area_ERA*Gv.ps_ERA[n,:,:] ) / np.sum( Gv.area_ERA )
        print( "ERA5 Global mean surface pressure=",globPS )
        
    #----------------------------------
    # Create a time array from on of
    # the ERA datasets
    #----------------------------------
    #pdTime_ERA = pd.to_datetime( dsT_ERA['time'].values )
    
    # Better/COnsistent to create time/date variables
    # For ERA5 let's add hour to ymdStr if hour0 != 99
    # Pandas understands space as delimiter between day
    # and hour
    if ( hour0 != 99 ):
        ymdhStr = ymdStr + ' ' + str( hour0 ).zfill(2)
    else:
        ymdhStr = ymdStr
    PdTime_ERA = pd.date_range( ymdhStr , periods=nt,freq='H')
    Gv.pdTime_ERA = pd.to_datetime( PdTime_ERA.values )
    

    #-----------------------------------------------
    # Get hybrid eta-coordinate coefficients for ERA5
    #-----------------------------------------------
    Gv.amid_ERA = dsT_ERA['a_model'].values 
    Gv.bmid_ERA = dsT_ERA['b_model'].values
    Gv.aint_ERA = dsT_ERA['a_half'].values 
    Gv.bint_ERA = dsT_ERA['b_half'].values
    print( "shape of a_model ", np.shape( Gv.amid_ERA ) ) 
    

    toc = time.perf_counter()
    pTime = f"Reading one set of ERA5 vars took  {toc - tic_overall:0.4f} seconds"
    print(pTime)

    # To make this code general there should be a split here,
    # i.e., after reading in all the data to process
    rcode=1
    return rcode


def get_ERAI( year=2022, month=11, day=1, hour0=99, interactive=False ):

    global pdTime_ERA
    global ps_ERA
    global te_ERA
    global q_ERA
    global u_ERA
    global v_ERA
    global w_ERA
    global amid_ERA, bmid_ERA, aint_ERA, bint_ERA
    # For diagnostic puroposes
    global lon_ERA, lat_ERA

    try:
        MySrc
        if ( MySrc != 'ERAI'):
            print( "You shouldnt be here - ABORT")
            rcode=-1
            return rcode
    except NameError:
        print( 'go on ' )

    tic_overall = time.perf_counter()

    # ERA style month string - yyyymm
    monStr=str( year ).zfill(4)+str(month).zfill(2)
    # CAM history style yyyy-mm-dd string
    ymdStr=str( year ).zfill(4) + '-' + str(month).zfill(2) + '-' + str(day).zfill(2)
    
    if ( hour0 != 99 ):
        hour1=hour0+5
        ymdh0=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+str(hour0).zfill(2)
        ymdh1=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+str(hour1).zfill(2)
        ymdh=ymdh0  
    else:
        ymdh0=str( year ).zfill(4)+str(month).zfill(2)+str(day).zfill(2)+"*"
        ymdh=ymdh0
        
    print( "Time tags for ERAI files ...")
    print(monStr)
    print(ymdh) 

    #wrkdir=  '/glade/scratch/juliob/erai_2017/' 
    wrkdir=  '/glade/scratch/juliob/erai_'+ str( year ).zfill(4) +'/' 

    #print(dstTZHkey)

    #Define all file names for later use in dask function
    #-----------------------------------------------------
    scfile = wrkdir + 'ei.oper.an.ml.regn128sc.'+ymdh+'.nc'
    uvfile = wrkdir + 'ei.oper.an.ml.regn128uv.'+ymdh+'.nc'
    
    print(scfile)
    print(uvfile)
    #Serial read
    #--------------
    dsSC_ERA   = xr.open_mfdataset( scfile, combine='nested',concat_dim=['time'], data_vars='different', coords='different' ) 
    dsUV_ERA   = xr.open_mfdataset( uvfile ,combine='nested',concat_dim=['time'], data_vars='different', coords='different' ) #data_vars='different', coords='different')

    ps_ERA = np.exp( dsSC_ERA['LNSP_GDS4_HYBL'].values )
    te_ERA = dsSC_ERA['T_GDS4_HYBL' ].values
    q_ERA  = dsSC_ERA['Q_GDS4_HYBL' ].values
    w_ERA  = dsSC_ERA['W_GDS4_HYBL' ].values
    u_ERA  = dsUV_ERA['U_GDS4_HYBL'].values
    v_ERA  = dsUV_ERA['V_GDS4_HYBL'].values
        
    lat_ERA = dsSC_ERA['g4_lat_0'].values
    lon_ERA = dsSC_ERA['g4_lon_1'].values

    print( " Checking ERA-i data ")
    print( " U-shape: ", np.shape(u_ERA) )
    print( " V-shape: ", np.shape(v_ERA) )
    print( " T-shape: ", np.shape(te_ERA) )
    print( " Q-shape: ", np.shape(q_ERA) )

    #---------------------------
    # Get shape of ERA data
    #-----------------------------
    nt,nL,nx,ny = np.shape( te_ERA )

    #----------------------------------
    # Create a time array from on of
    # the ERA datasets
    #----------------------------------
    
    PdTime_ERA = pd.date_range( ymdStr , periods=nt,freq='6H')
    pdTime_ERA = pd.to_datetime( PdTime_ERA.values )
    #-----------------------------------------------
    # Get hybrid eta-coordinate coefficients for ERAI
    # Why do they have different names for
    # these in the 'uv' and 'sc' files????
    #-----------------------------------------------
    amid_ERA = dsSC_ERA['lv_HYBL2_a'].values 
    bmid_ERA = dsSC_ERA['lv_HYBL2_b'].values
    aint_ERA = dsSC_ERA['lv_HYBL_i3_a'].values 
    bint_ERA = dsSC_ERA['lv_HYBL_i3_b'].values
    print( "shape of hybrid a in ERAI ", np.shape( amid_ERA ) ) 




    toc = time.perf_counter()
    pTime = f"Reading one set of ERA5 vars took  {toc - tic_overall:0.4f} seconds"
    print(pTime)

    # To make this code general there should be a split here,
    # i.e., after reading in all the data to process
    rcode=1

    return rcode

def get_Src( year=2022, month=11, day=1, hour0=99):
    
    if (Gv.MySrc=='ERA5'):
        faa =get_ERA5( year=year, month=month, day=day, hour0=hour0)
        rcode=faa
    else:
        print('Aaaaaaaaaaaa')
        rcode=-999
        
    return rcode
