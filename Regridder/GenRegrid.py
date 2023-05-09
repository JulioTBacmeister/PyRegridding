#!/usr/bin/env python
# Import packages 
import sys
import argparse as arg

import xarray as xr
import numpy as np
import pandas as pd

from scipy.io import FortranFile
from scipy import interpolate as intr
from scipy.interpolate import RegularGridInterpolator as RGi


import ESMF as E

import importlib
import glob
import copy
import time

import scripGen as SG
import esmfRegrid as erg

#import cfgrib


import dask
import dask.array as da

#import VertRegrid as vrg

# "ChatGPI version" --- 
import VertRegridFlexLL as vrg
print( "Using Flexible parallel/serial VertRegrid ")

import FVRegrid as FV

# import modules in other directories
sys.path.append('../Utils/')
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

def prep(Dst = 'ne30pg3', DstVgrid='L58',  Src='ERA5', WOsrf=False , RegridMethod="CONSERVE" ):
    #---------------------------------------------
    # This function sets-up variables and objects 
    # that are need for horizontal and vertical 
    # regridding of ERA reanalyses.
    #---------------------------------------------

    global MyDst, MyDstVgrid, MySrc
    global regrd,srcf,dstf
    global phis_CAM, phis_ERA
    global amid_CAM,bmid_CAM,aint_CAM,bint_CAM
    global lon_CAM,lat_CAM,area_CAM
    global area_ERA

    # Grid keys for remapping
    global srcHkey,dstHkey,srcTHkey,dstTHkey,srcZHkey,dstZHkey,srcTZHkey,dstTZHkey

    global doWilliamsonOlson

    global p_00_ERA, p_00_CAM

    #------- 
    # Begin
    #-------
    tic_overall = time.perf_counter()
    MyDst,MyDstVgrid,MySrc = Dst,DstVgrid,Src

    doWilliamsonOlson = WOsrf
    p_00_CAM = 100_000.

    print( f"In prep Src= {Src} to Dst={Dst} " )
    if (Dst == 'ne30pg3'):
        dstHkey = 'c'
        dst_type='mesh'
        dst_scrip = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/ne30pg3_scrip_170611.nc'
        dst_TopoFile = '/glade/p/cgd/amp/juliob/bndtopo/latest/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_20230105.nc'

    if ((Dst == 'fv0.9x1.25') or (Dst=='fv1x1')):
        dstHkey = 'yx'
        dst_type='grid'
        dst_scrip = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.9x1.25_141008.nc'
        #dst_TopoFile='/glade/p/cgd/amp/juliob/bndtopo/latest/fv_0.9x1.25_gmted2010_modis_bedmachine_nc3000_Laplace0100_20220708.nc'
        dst_TopoFile = '/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_160505.nc'

    if (Src == 'ERA5'):
        srcHkey = 'yx'
        src_type='grid'
        src_scrip = '/glade/work/juliob/ERA5-proc/ERA5interp/grids/ERA5_640x1280_scrip.nc'
        src_TopoFile = '/glade/work/juliob/ERA5-proc/ERA5interp/phis/ERA5_phis.nc'
        p_00_ERA = 1.0

    if (Src == 'ERAI'):
        srcHkey = 'yx'
        src_type='grid'
        #print(" READING PRECOMPUTED SCRIP FROM JERRY/PATRICK CODE !!!!!!!! " )
        #src_scrip = '/glade/work/juliob/ERA-I-grids/Sgrid_SRC.nc'
        src_scrip = '/glade/work/juliob/ERA-I-grids/ERAI_256x512_scrip.nc'
        src_TopoFile = '/glade/scratch/juliob/erai_2017/ei.oper.an.ml.regn128sc.2017010100.nc'
        p_00_ERA = 100_000.

    # ----------------------------------------------
    # Get DST vertical grid from a file.
    # These should be small files but 
    # they aren't always.
    # ----------------------------------------------
    if (DstVgrid == 'L58' ):
        # Read in CAM L58 vertical grid
        dstVgridFile = '/glade/work/juliob/ERA5-proc/CAM-grids/Vertical/GRID_48_taperstart10km_lowtop_BL10_v3p1_beta1p75.nc'

    if (DstVgrid == 'L32' ):
        dstVgridFile = '/glade/p/cesmdata/cseg/inputdata/atm/cam/inic/se/f.e22.FC2010climo.ne30pg3_ne30pg3_mg17.cam6_2_022.002.cam.i.0020-01-01-00000_c200610.nc'

    # Set grid keys for Src ERA5 reanalysis
    srcTHkey  = 't'  + srcHkey
    srcZHkey  = 'z'  + srcHkey
    srcTZHkey = 'tz' + srcHkey

    # Set grid keys for Dst CAM-SE
    dstTHkey  = 't'  + dstHkey
    dstZHkey  = 'z'  + dstHkey
    dstTZHkey = 'tz' + dstHkey
 
    # ----------------------------------------------
    # Get all topo data we will use
    # Read in CAM topography. Also get
    # lon and lat and area for CAM (Dst)
    # grid.
    # ----------------------------------------------
    dsTopo_CAM=xr.open_dataset( dst_TopoFile )
    varsCAM  = list( dsTopo_CAM.variables )
    phis_CAM = dsTopo_CAM['PHIS'].values
    lon_CAM  = dsTopo_CAM['lon'].values
    lat_CAM  = dsTopo_CAM['lat'].values
    if ('area' in varsCAM):
        area_CAM = dsTopo_CAM['area'].values
    else:
        area_CAM = GrU.area2d( lon=lon_CAM, lat=lat_CAM )

    if (Src == 'ERA5'):
        # Read in ERA5 topography
        dsTopo_ERA=xr.open_dataset( src_TopoFile )
        phis_ERA=dsTopo_ERA['Z_GDS4_SFC'].values

    if (Src == 'ERAI'):
        # Read in ERA-I topography
        dsTopo_ERA=xr.open_dataset( src_TopoFile )
        phis_ERA=dsTopo_ERA['Z_GDS4_HYBL'].values

    # ----------------------------------------------
    # Look for pre-computed weights file
    # If none, set params to create weights file
    # ----------------------------------------------
    if ( (Src == 'ERA5') and (Dst == 'ne30pg3') ):
        griddir = "/glade/work/juliob/ERA5-proc/ERA5interp/grids/"
        wgts_file_Con = griddir + "ERA5_ne30pg3_Conserv_wgts.nc"
        write_weights = False 
        read_weights = True 
    else:
        wgts_file_Con = "REGRID_"+Src+"_x_"+Dst+"_"+RegridMethod+".nc"
        write_weights = False 
        read_weights = False 


    """
    elif ( (Src == 'ERAI') and (Dst == 'fv0.9x1.25') or (Dst=='fv1x1') ):
        print(" READING PRECOMPUTED WEIGHTS FROM JERRY/PATRICK CODE !!!!!!!! " )
        griddir = "/glade/work/juliob/ERA-I-grids/"
        wgts_file_Con = griddir + "Swgt_SRC2DST.nc"
        write_weights = False 
        read_weights = True 
    """


    # ----------------------------------------------
    #  Set-up regridding machinery
    # ----------------------------------------------
    # Scrip file for ERA5 created by ERA5scrip.ipynb
    if (Src == 'ERA5'):
        dsERAscrip = xr.open_dataset( src_scrip )
        area_ERA = -9999. #np.reshape( dsERAscrip['grid_area'].values , np.shape( phis_ERA ) )
    else:
        area_ERA = -9999.
        
    # ----------------------------------------------
    # Make object for ESMF regridding from SRC
    # grid to CAM target. Scrip files need to be provided even 
    # when a weight file is used
    # ----------------------------------------------
    regrd, srcf, dstf = erg.Regrid( srcScrip = src_scrip , 
                                    srcType  = src_type  ,
                                    dstScrip = dst_scrip ,
                                    dstType  = dst_type  ,
                                    write_weights = write_weights ,
                                    read_weights = read_weights ,
                                    weights_file = wgts_file_Con ,
                                    RegridMethod = RegridMethod )
    


    vCAM=xr.open_dataset( dstVgridFile )
    amid_CAM = vCAM['hyam'].values
    bmid_CAM = vCAM['hybm'].values
    aint_CAM = vCAM['hyai'].values
    bint_CAM = vCAM['hybi'].values

    print( f" Src scripfile {src_scrip} " )
    print( f" Dst scripfile {dst_scrip} " )
    print( f" Src topo file {src_TopoFile} " )
    print( f" Dst topo file {dst_TopoFile} " )
    print( f" {DstVgrid} Dst vertical grid from {dstVgridFile} " )


    toc = time.perf_counter()
    pTime = f"Prepping for {Src} to {Dst} proc in {__name__} took  {toc - tic_overall:0.4f} seconds"
    print(pTime)
 
    code = 1
    return code


# Define a function that loads a NetCDF file and returns an xarray dataset
@dask.delayed
def load_file(path):
    ds = xr.open_mfdataset(path ,  data_vars='different', coords='different' )
    return ds


def get_ERA5( year=2022, month=11, day=1, hour0=99):

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
        if ( MySrc != 'ERA5'):
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

    era5dir = "/glade/collections/rda/data/ds633.6/e5.oper.an.ml/"
    wrkdir=era5dir+monStr+"/"

    print(dstTZHkey)

    #Define all file names for later use in dask function
    #-----------------------------------------------------
    spfile= wrkdir + 'e5.oper.an.ml.128_134_sp.regn320sc.'+ymdh+'.nc'
    tfile = wrkdir + 'e5.oper.an.ml.0_5_0_0_0_t.regn320sc.'+ymdh+'.nc'
    qfile = wrkdir + 'e5.oper.an.ml.0_5_0_1_0_q.regn320sc.'+ymdh+'.nc'
    ufile = wrkdir + 'e5.oper.an.ml.0_5_0_2_2_u.regn320uv.'+ymdh+'.nc'
    vfile = wrkdir + 'e5.oper.an.ml.0_5_0_2_3_v.regn320uv.'+ymdh+'.nc'
    wfile = wrkdir + 'e5.oper.an.ml.0_5_0_2_8_w.regn320sc.'+ymdh+'.nc'
    all_ERA_files = [ spfile , tfile, qfile, ufile, vfile, wfile ]

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
    dsPS_ERA  = xr.open_mfdataset( spfile, data_vars='different', coords='different' )
    dsT_ERA   = xr.open_mfdataset( tfile , data_vars='different', coords='different')
    dsQ_ERA   = xr.open_mfdataset( qfile,  data_vars='different', coords='different' )
    dsU_ERA   = xr.open_mfdataset( ufile , data_vars='different', coords='different')
    dsV_ERA   = xr.open_mfdataset( vfile , data_vars='different', coords='different')
    dsW_ERA   = xr.open_mfdataset( wfile , data_vars='different', coords='different')
    """

    ps_ERA = dsPS_ERA['SP'].values
    te_ERA = dsT_ERA['T'].values
    q_ERA  = dsQ_ERA['Q'].values
    u_ERA  = dsU_ERA['U'].values
    v_ERA  = dsV_ERA['V'].values
    w_ERA  = dsW_ERA['W'].values

    lon_ERA = dsT_ERA['longitude'].values
    lat_ERA = dsT_ERA['latitude'].values
    
    #---------------------------
    # Get shape of ERA data
    #-----------------------------
    nt,nL,nx,ny = np.shape( te_ERA )

    #    Have a look at ERA global mean surface pressures
    for n in np.arange(nt):
        globPS = np.sum( area_ERA*ps_ERA[n,:,:] ) / np.sum( area_ERA )
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
    pdTime_ERA = pd.to_datetime( PdTime_ERA.values )
    

    #-----------------------------------------------
    # Get hybrid eta-coordinate coefficients for ERA5
    #-----------------------------------------------
    amid_ERA = dsT_ERA['a_model'].values 
    bmid_ERA = dsT_ERA['b_model'].values
    aint_ERA = dsT_ERA['a_half'].values 
    bint_ERA = dsT_ERA['b_half'].values
    print( "shape of a_model ", np.shape( amid_ERA ) ) 
    

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

    """
    #    Have a look at ERA global mean surface pressures
    #    Need to get an area for ERA-I grid.
    for n in np.arange(nt):
        globPS = np.sum( area_ERA*ps_ERA[n,:,:] ) / np.sum( area_ERA )
        print( "ERA5 Global mean surface pressure=",globPS )
    """

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

    if (interactive==True):
        print( "Interactive !!!! " )
        print( pdTime_ERA )
        return pdTime_ERA, dsSC_ERA
    else:
        return rcode

def xRegrid( ExitAfterTemperature=False , 
             ExitAfterWinds=False , 
             HorzInterpLnPs=False , 
             Use_ps_ERA_xCAM_in_vert=True ):

    # Variable created in this function and needed elsehwere in
    # the module
    #----------------------------------------------------------
    global ps_CAM          
    global phis_ERA_xCAM
    global te_ERA_xzCAM
    global q_ERA_xzCAM
    global u_ERA_xzCAM
    global v_ERA_xzCAM
    global w_ERA_xzCAM

    StartTime = time.asctime( time.localtime(time.time()) )
    tic_overall = time.perf_counter()
    print( f"starting xRegrid {MySrc} _x_ {MyDst} at {StartTime} ")
    #-----------------------------------------
    # Horz remap of PHIS - No time dimension
    #-----------------------------------------
    phis_ERA_xCAM = erg.HorzRG( aSrc = phis_ERA , 
                                regrd = regrd , 
                                srcField=srcf , 
                                dstField=dstf , 
                                srcGridkey=srcHkey ,
                                dstGridkey=dstHkey )
    
    toc_here = time.perf_counter()
    pTime = f"Finished phis Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)

    """
    #------------------------------------------
    # Step in with Islas's 1440x720 =>xCAM 
    # f09 regridded ERA5 topo.
    # This is for diagnostics/debugging only
    #-----------------------------------------
    fTopo2='/glade/u/home/islas/for/suqin/regridera5/makephis/output/PHIS_model_and_ERA5_analysis_f09_f09.nc'
    dTopo2=xr.open_dataset( fTopo2 )
    phis_ERA_xCAM =dTopo2['PHIS_analysis'].values[0,:,:]
    """
 
    
    #-----------------------------------------
    # Calculate difference between phis's
    #-----------------------------------------
    Dphis = phis_ERA_xCAM - phis_CAM
    
        
    #-----------------------------------------
    # Horz remap of PS
    # log-exp bit is to reproduce W&O
    #-----------------------------------------
    if (HorzInterpLnPs==True):
        xfld_ERA = np.log( ps_ERA )
    else:
        xfld_ERA = ps_ERA
        
    ps_ERA_xCAM    = erg.HorzRG( aSrc = xfld_ERA , 
                                 regrd = regrd , 
                                 srcField=srcf , 
                                 dstField=dstf , 
                                 srcGridkey=srcTHkey ,
                                 dstGridkey=dstHkey )

    if (HorzInterpLnPs==True):
        ps_ERA_xCAM = np.exp( ps_ERA_xCAM ) 

    toc_here = time.perf_counter()
    pTime = f"Finished ps Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)

    #-----------------------------------------
    # Make 3D prssure field on ERA ZH grid
    #-----------------------------------------
    pmid_ERA, pint_ERA, delp_ERA \
        = MkP.Pressure (am=amid_ERA ,
                        bm=bmid_ERA ,
                        ai=aint_ERA ,
                        bi=bint_ERA ,
                        ps=ps_ERA ,
                        p_00=p_00_ERA ,
                        Gridkey = srcTZHkey )

    #-----------------------------------------
    # Horz remap of temperature
    #-----------------------------------------
    te_ERA_xCAM    = erg.HorzRG( aSrc = te_ERA , 
                                 regrd = regrd , 
                                 srcField=srcf , 
                                 dstField=dstf , 
                                 srcGridkey=srcTZHkey ,
                                 dstGridkey= dstHkey )
    
    toc_here = time.perf_counter()
    pTime = f"Finished te_ERA Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)

    #-----------------------------------------------------------
    # Find "T_bot" and "P_bot" defined in Williamson & Olson 
    # as occurring at "the first level above 150m ... " 
    # This uses ERA temperature and surface pressure remapped  
    # to CAM horz grid (adapted from Isla's function)
    #-----------------------------------------------------------
    


    print( "amid_ERA    ", np.shape( amid_ERA ) )
    print( "ps_ERA_xCAM ", np.shape( ps_ERA_xCAM ) )
    print( "te_ERA_xCAM ", np.shape( te_ERA_xCAM ) )

    tic_150 = time.perf_counter()

    te_150, pmid_150,L150  = MkP.Pressure_TandP150 ( am=amid_ERA ,
                                                     bm=bmid_ERA ,
                                                     ai=aint_ERA ,
                                                     bi=bint_ERA ,
                                                     ps=ps_ERA_xCAM ,
                                                     te=te_ERA_xCAM , 
                                                     p_00=p_00_ERA , 
                                                     Gridkey = dstTZHkey )


    toc_150 = time.perf_counter()
    pTime = f"Finding Te150 and P150 took  {toc_150 - tic_150:0.4f} seconds"
    print(pTime)

    
    #-------------------------------------------------------------------
    #                    "CAM surface pressure"
    #-------------------------------------------------------------------
    # We don't actually have ps from CAM, so we make a guess based on
    # ps_ERA_xCAM and te_bot asdescribed above.  In a sense this is a 
    # vertical remapping, so we could call this variable ps_ERA_xzCAM, 
    # but we'll just call it ps_CAM ...
    #-------------------------------------------------------------------
    
    
    ps_CAM = vrg.PsAdjust( phis=phis_ERA_xCAM, 
                           phis_CAM=phis_CAM, 
                           ps=ps_ERA_xCAM , 
                           pm150=pmid_150 , 
                           te150=te_150 , 
                           Gridkey=dstTZHkey  )

    """
    
    gravit  = 9.80616           # acceleration of gravity ~ m/s^2
    boltz   = 1.38065e-23       # boltzmann's constant ~ J/k/molecule
    avogad  = 6.02214e26        # avogadro's number ~ molecules/kmole
    mwdair  = 28.966            # molecular weight dry air ~ kg/kmole

    rgas    = avogad*boltz      # universal gas constant ~ J/k/kmole
    rdair   = rgas/mwdair       #constant for dry air   ~ J/k/kg

    if (dstTZHkey == 'tzc'):
        nt,nz,ncol = np.shape(te_ERA_xCAM)
        L_bot = nz-1
        te_bot = te_150 # te_ERA_xCAM[:,L_bot,:]
        ps_CAM = np.zeros( (nt , ncol) )
        for i in np.arange( nt ):
            ps_CAM[i,:] = ps_ERA_xCAM[i,:] * np.exp( Dphis / (rdair*te_bot[i,:]) )

    if (dstTZHkey == 'tzyx'):
        nt,nz,ny,nx = np.shape(te_ERA_xCAM)
        L_bot = nz-1
        #te_bot = te_150 # te_ERA_xCAM[:,L_bot,:]
        print( "using simpl TBOT fro ps_adj !!!!!! ")
        te_bot = te_ERA_xCAM[:,L_bot,:]
        ps_CAM = np.zeros( (nt , ny, nx) )
        for i in np.arange( nt ):
            ps_CAM[i,:,:] = ps_ERA_xCAM[i,:,:] * np.exp( Dphis / (Rdry*te_bot[i,:,:]) )
    
    
    print( " RDAIR in w and O", rdair , " NZ here ",nz )
    print( " RDry here ", Rdry  )
    """
    
    #-----------------------------------------------------------------------------------------------
    # Now we creat full 4(3)D pressure fields on the ERA and CAM vertical grids. These are used for 
    # vertical interpolation below. Not clear what surface pressure to use when building ERA vertical
    # grid, i.e., ps_CAM or ps_ERA_xCAM.
    #
    # The WO2015 document seems to suggest ps_CAM but W&O code definitely uses ps_ERA_xCAM
    #-----------------------------------------------------------------------------------------------
    tic_P3D = time.perf_counter()
    
    # Surface pressure: Choose wisely
    #---------------------------------
    if (Use_ps_ERA_xCAM_in_vert == True):
        ps_FLD = ps_ERA_xCAM
    else:
        ps_FLD = ps_CAM
        
    pmid_CAM_zERA, pint_CAM_zERA, delp_CAM_zERA \
        = MkP.Pressure (am=amid_ERA ,
                        bm=bmid_ERA ,
                        ai=aint_ERA ,
                        bi=bint_ERA ,
                        ps=ps_FLD , # What to use here seems key: ps_CAM or ps_ERA_xCAM
                        p_00=p_00_ERA ,
                        Gridkey = dstTZHkey )

    pmid_CAM,pint_CAM,delp_CAM \
        = MkP.Pressure (am=amid_CAM ,
                        bm=bmid_CAM ,
                        ai=aint_CAM ,
                        bi=bint_CAM ,
                        ps=ps_CAM ,
                        p_00=p_00_CAM , 
                        Gridkey = dstTZHkey )

    #-----------------------------------------------------------
    # Log-pressure is preferred for vertical interpolation
    # per Williamson&Olson
    #-----------------------------------------------------------
    p_00 = 100_000. # Here we just use the sensible value of p_00
    lnpint_CAM = -7_000. * np.log( pint_CAM / p_00 )
    lnpmid_CAM = -7_000. * np.log( pmid_CAM / p_00 )
    lnpmid_CAM_zERA = -7_000. * np.log( pmid_CAM_zERA /p_00 )

    toc_P3D = time.perf_counter()
    pTime = f"Creating 3D P-fields etc., took   {toc_P3D - tic_P3D:0.4f} seconds"
    print(pTime)
    
    
    if ( doWilliamsonOlson == True ):
        tic_WO = time.perf_counter()
        print( "WilliamsonOlson surface " )
        #----------------------------------------------------------
        # Calculate extrapolated surface temperature using 
        # Williamson & Olson standard lapse rate approach
        #----------------------------------------------------------
        ts_extrap = vrg.TsExtrap( ps = ps_CAM ,
                                  pm150 = pmid_150 ,
                                  te150 = te_150 )

                              
        #--------------------------------------------------------
        # If Williamson & Olson treatment of surface layer 
        # is selected then correct te_ERA_xzCAM between 
        # pmid_150 and ps_CAM
        #-------------------------------------------------------
        te_WO  =    vrg.TeWO( te = te_ERA_xCAM ,
                              pmid = pmid_CAM_zERA ,
                              te150 = te_150 ,
                              pm150 = pmid_150 ,
                              ts = ts_extrap,
                              ps = ps_CAM, 
                              L150 = L150 ,
                              Gridkey =dstTZHkey )

        te_ERA_xCAM = copy.deepcopy( te_WO )
        toc_WO = time.perf_counter()
        pTime = f"Williamson Olson surface took  {toc_WO - tic_WO:0.4f} seconds"
        print(pTime)                       
    

    print(" going into vertical regrid of T " )
    te_ERA_xzCAM = vrg.VertRG( a_x  = te_ERA_xCAM ,
                               zSrc = lnpmid_CAM_zERA ,
                               zDst = lnpmid_CAM ,
                               Gridkey =dstTZHkey ,
                               kind = 'quadratic' ) #linea
    

    #-------------------------------------------------------------
    # return statement to assist in debugging and analysis
    #-------------------------------------------------------------
    if ( ExitAfterTemperature == True ):
        return pmid_ERA, lat_ERA, lon_ERA, te_ERA, \
            pmid_CAM_zERA, lat_CAM, lon_CAM, te_ERA_xCAM, \
            pmid_CAM, te_ERA_xzCAM, \
            ps_ERA, ps_CAM, ps_ERA_xCAM, \
            phis_ERA, phis_CAM, phis_ERA_xCAM



    #-------------------------------------------------------------
    # Continue on to regridding Q, U, V, W ...
    #-------------------------------------------------------------

    #--------------------
    #  Regridding of Q
    #---------------------
    print(" going into horz+vertical regrid of Q " )
    q_ERA_xzCAM , q_ERA_xCAM = fullRegrid ( a_ERA = q_ERA ,
                                            zSrc = lnpmid_CAM_zERA ,
                                            zDst = lnpmid_CAM )
        
    q_ERA_xzCAM = vrg.BottomFill( a_zCAM = q_ERA_xzCAM ,
                                  a_zERA = q_ERA_xCAM ,
                                  pmid_zCAM=pmid_CAM ,
                                  ps_ERA = ps_ERA_xCAM , 
                                  Gridkey = dstTZHkey )

    qx = SaturateQ( q=q_ERA_xzCAM , 
                    te=te_ERA_xzCAM ,
                    p=pmid_CAM, 
                    Gridkey = dstTZHkey )
    
    q_ERA_xzCAM =  copy.deepcopy(qx)




    #--------------------
    #  Regridding of U
    #---------------------
    print(" going into horz+vertical regrid of U " )
    u_ERA_xzCAM, u_ERA_xCAM = fullRegrid ( a_ERA = u_ERA ,
                                           zSrc = lnpmid_CAM_zERA ,
                                           zDst = lnpmid_CAM )

    u_ERA_xzCAM = vrg.BottomFill( a_zCAM = u_ERA_xzCAM ,
                                  a_zERA = u_ERA_xCAM ,
                                  pmid_zCAM=pmid_CAM ,
                                  ps_ERA = ps_ERA_xCAM , 
                                  Gridkey = dstTZHkey )
    
    #--------------------
    #  Regridding of V
    #---------------------
    print(" going into horz+vertical regrid of V " )
    v_ERA_xzCAM, v_ERA_xCAM = fullRegrid ( a_ERA = v_ERA ,
                                           zSrc = lnpmid_CAM_zERA ,
                                           zDst = lnpmid_CAM )

    v_ERA_xzCAM = vrg.BottomFill( a_zCAM = v_ERA_xzCAM ,
                                  a_zERA = v_ERA_xCAM ,
                                  pmid_zCAM=pmid_CAM ,
                                  ps_ERA = ps_ERA_xCAM , 
                                  Gridkey = dstTZHkey )
    
    #--------------------
    #  Regridding of W
    #---------------------
    print(" going into horz+vertical regrid of W " )
    w_ERA_xzCAM, w_ERA_xCAM = fullRegrid ( a_ERA = w_ERA ,
                                          zSrc = lnpmid_CAM_zERA ,
                                          zDst = lnpmid_CAM )

    w_ERA_xzCAM = vrg.BottomFill( a_zCAM = w_ERA_xzCAM ,
                                  a_zERA = w_ERA_xCAM ,
                                  pmid_zCAM=pmid_CAM ,
                                  ps_ERA = ps_ERA_xCAM , 
                                  Gridkey = dstTZHkey )
    





    #-------------------------------------------------------------
    # return statement to assist in debugging and analysis
    #-------------------------------------------------------------
    if ( ExitAfterWinds == True ):
        return pmid_ERA, lat_ERA, lon_ERA, u_ERA, \
            pmid_CAM_zERA, lat_CAM, lon_CAM, u_ERA_xCAM, \
            pmid_CAM, u_ERA_xzCAM, \
            ps_ERA, ps_CAM, ps_ERA_xCAM, \
            phis_ERA, phis_CAM, phis_ERA_xCAM



    
    toc = time.perf_counter()
    pTime = f"Overall time in this function  {toc - tic_overall:0.4f} seconds"
    print(pTime)
        
    rcode =1 
    return rcode

def fullRegrid( a_ERA,  zSrc ,  zDst , kind='linear', ReturnVars=2 ):
    
    print("Horz RG in fullRegrid " )
    a_ERA_xCAM    = erg.HorzRG( aSrc = a_ERA , 
                                 regrd = regrd , 
                                 srcField=srcf , 
                                 dstField=dstf , 
                                 srcGridkey=srcTZHkey,
                                 dstGridkey=dstHkey )

    print("Vert RG in fullRegrid " )
    a_ERA_xzCAM = vrg.VertRG( a_x  = a_ERA_xCAM ,
                              zSrc = zSrc ,
                              zDst = zDst ,
                              Gridkey=dstTZHkey,
                              kind = kind )

    if (ReturnVars==1):
        return a_ERA_xzCAM
    if (ReturnVars==2):
        return a_ERA_xzCAM,a_ERA_xCAM

def SaturateQ ( q , te , p, Gridkey ):

    qsat = hum.qsat( p=p, T=te )
    qx=np.minimum( q , qsat )

    """
    print( "UGLY temporary Kluge: Don't saturate/ZERO top 3 levels of model " )
    if (Gridkey == 'tzyx' ):
        qx[:,0:2,:,:]= 0*q[:,0:2,:,:]
    if (Gridkey == 'tzc' ):
        qx[:,0:2,:]= 0*q[:,0:2,:]
    if (Gridkey == 'zyx' ):
        qx[0:2,:,:]= 0*q[0:2,:,:]
    if (Gridkey == 'zc' ):
        qx[0:2,:]= 0*q[0:2,:]
    
    """
    
    return qx

def write_netcdf( version='' ):
    ntime = np.shape(pdTime_ERA)[0]
    print(ntime)
    Bfilo="/glade/scratch/juliob/" + MySrc +"_x_"+ MyDst + "_"+ MyDstVgrid + "_" + version #  + "."     #  + timetag+ ".nc"

    if (doWilliamsonOlson==True):
        Bfilo = Bfilo + '_WO'

    ilev = (aint_CAM+bint_CAM)* 1_000. #* 100_000.
    lev  = (amid_CAM+bmid_CAM)* 1_000. #* 100_000.



    if (dstTZHkey == 'tzc' ):
        nt,nz,ncol = np.shape( te_ERA_xzCAM )
        for itim in np.arange( ntime ):
            dims   = ["ncol","time","lev","ilev"]
            coords = dict( 
                lon  = ( ["ncol"],lon_CAM ),
                lat  = ( ["ncol"],lat_CAM ),
                lev  = ( ["lev"],lev),
                ilev = ( ["ilev"],ilev),
                time = ( ["time"],  np.array(itim ,ndmin=1 ) ), #pd.to_datetime( pdTime_ERA[itim] ) ),
            )
        
            Wds = xr.Dataset( coords=coords  )
            Wds["TimeStamp"] = pd.to_datetime( pdTime_ERA[itim] )
            Wds["P_00"] = 100_000.
        
            Dar = xr.DataArray( data=aint_CAM, dims=('ilev',),
                                attrs=dict( description='interface hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyai'] = Dar

            Dar = xr.DataArray( data=bint_CAM, dims=('ilev',),
                                attrs=dict( description='interface hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybi'] = Dar

            Dar = xr.DataArray( data=amid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyam'] = Dar

            Dar = xr.DataArray( data=bmid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybm'] = Dar
        
            Dar = xr.DataArray( data=area_CAM, dims=('ncol',),
                                attrs=dict( description='Cell area',units='Steradians',) ,) 
            Wds['area'] = Dar

            Dar = xr.DataArray( data=phis_CAM, dims=('ncol',),
                                attrs=dict( description='Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS'] = Dar

            Dar = xr.DataArray( data=phis_ERA_xCAM, dims=('ncol',),
                                attrs=dict( description='ERA Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS_ERA'] = Dar

            Dar = xr.DataArray( data=ps_CAM[itim,:].reshape(1,ncol), 
                                dims=('time','ncol',),
                                attrs=dict( description='Surface Pressure',units='Pa',) ,) 
            Wds['PS'] = Dar
    
            Dar = xr.DataArray( data=te_ERA_xzCAM[itim,:,:].reshape(1,nz,ncol), 
                                dims=('time','lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['T'] = Dar

            Dar = xr.DataArray( data=q_ERA_xzCAM[itim,:,:].reshape(1,nz,ncol), 
                                dims=('time','lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['Q'] = Dar
        
            Dar = xr.DataArray( data=u_ERA_xzCAM[itim,:,:].reshape(1,nz,ncol), 
                                dims=('time','lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['U'] = Dar

            Dar = xr.DataArray( data=v_ERA_xzCAM[itim,:,:].reshape(1,nz,ncol), 
                                dims=('time','lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['V'] = Dar

            Dar = xr.DataArray( data=w_ERA_xzCAM[itim,:,:].reshape(1,nz,ncol), 
                                dims=('time','lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['W'] = Dar

            yymmdd = str(pdTime_ERA[itim])[0:10]
            hr=str(pdTime_ERA[itim])[11:13]
            ss = str(int(hr)*3600).zfill(5)
            timetag =  yymmdd+'-'+ss
            filo= Bfilo + "." + timetag+ ".nc"
            print( filo )
            Wds.to_netcdf( filo ,format="NETCDF3_CLASSIC" )

    if (dstTZHkey == 'tzyx' ):
        tic_FVstag = time.perf_counter()
        US,VS,slat,slon = FV.uvStaggers(U=u_ERA_xzCAM, 
                                        V=v_ERA_xzCAM,
                                        lon=lon_CAM,
                                        lat=lat_CAM   )
        toc_FVstag = time.perf_counter()
        pTime = f"Creating FV staggered US,VS took {toc_FVstag - tic_FVstag:0.4f} seconds"
        print(pTime)
        
        nt,nz,ny,nx = np.shape( te_ERA_xzCAM )
        for itim in np.arange( ntime ):
            dims   = ["lon","lat","time","lev","ilev"]
            coords = dict( 
                lon  = ( ["lon"],lon_CAM ),
                lat  = ( ["lat"],lat_CAM ),
                slon  = ( ["slon"],slon ),
                slat  = ( ["slat"],slat ),
                lev  = ( ["lev"],lev),
                ilev = ( ["ilev"],ilev),
                time = ( ["time"],  np.array(itim ,ndmin=1 ) ), #pd.to_datetime( pdTime_ERA[itim] ) ),
            )
        
            Wds = xr.Dataset( coords=coords  )
            Wds["TimeStamp"] = pd.to_datetime( pdTime_ERA[itim] )
            Wds["P_00"] = 100_000.
        
            Dar = xr.DataArray( data=aint_CAM, dims=('ilev',),
                                attrs=dict( description='interface hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyai'] = Dar

            Dar = xr.DataArray( data=bint_CAM, dims=('ilev',),
                                attrs=dict( description='interface hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybi'] = Dar

            Dar = xr.DataArray( data=amid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyam'] = Dar

            Dar = xr.DataArray( data=bmid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybm'] = Dar
        
            Dar = xr.DataArray( data=area_CAM, dims=('lat','lon',),
                                attrs=dict( description='Cell area',units='Steradians',) ,) 
            Wds['area'] = Dar

            Dar = xr.DataArray( data=phis_CAM, dims=('lat','lon',),
                                attrs=dict( description='Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS'] = Dar

            Dar = xr.DataArray( data=phis_ERA_xCAM, dims=('lat','lon',),
                                attrs=dict( description='ERA Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS_ERA'] = Dar

            Dar = xr.DataArray( data=ps_CAM[itim,:,:].reshape(1,ny,nx) , 
                                dims=('time','lat','lon',),
                                attrs=dict( description='Surface Pressure',units='Pa',) ,) 
            Wds['PS'] = Dar
    
            Dar = xr.DataArray( data=te_ERA_xzCAM[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['T'] = Dar

            Dar = xr.DataArray( data=q_ERA_xzCAM[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['Q'] = Dar
        
            Dar = xr.DataArray( data=u_ERA_xzCAM[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['U'] = Dar

            Dar = xr.DataArray( data=v_ERA_xzCAM[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['V'] = Dar

            Dar = xr.DataArray( data=US[itim,:,:,:].reshape(1,nz,ny-1,nx), 
                                dims=('time','lev','slat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['US'] = Dar

            Dar = xr.DataArray( data=VS[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','slon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['VS'] = Dar

            Dar = xr.DataArray( data=w_ERA_xzCAM[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['W'] = Dar
        
            yymmdd = str(pdTime_ERA[itim])[0:10]
            hr=str(pdTime_ERA[itim])[11:13]
            ss = str(int(hr)*3600).zfill(5)
            timetag =  yymmdd+'-'+ss
            filo= Bfilo + "." + timetag+ ".nc"
            print( filo )
            Wds.to_netcdf( filo ,format="NETCDF3_CLASSIC" )

    code = 1
    return code

#def Driver():
def main(year,month,day,hour):
    import calendar
    
    tic_total = time.perf_counter()
    days_in_month = calendar.monthrange(year,month)[1]

    print( f"About to process {year:n}-{month:n}-{day:n}")



    RegridMethod = 'CONSERVE'

    """
    DstVgrid='L58'
    Dst='ne30pg3'
    Src='ERA5'
    """
    DstVgrid='L32'
    Dst='fv1x1'
    Src='ERA5'

    lnPS=False
    if(lnPS==True):
        ver='lnPS'
    else:
        ver=''


    ret1 = prep(Dst=Dst, DstVgrid=DstVgrid ,Src=Src, WOsrf=True, RegridMethod=RegridMethod )
    if (day==99):
        for iday in np.arange( days_in_month):
            ret2 = get_ERA5( year=year ,month=month ,day=iday+1 , hour0=99 )
            ret3 = xRegrid(HorzInterpLnPs=lnPS )
            ret4 = write_netcdf(version=ver+'Test01')
    else:
        ret2 = get_ERA5( year=year ,month=month ,day=day , hour0=hour )
        ret3 = xRegrid(HorzInterpLnPs=lnPS )
        ret4 = write_netcdf(version=ver+'Test01')
        
    code = 1
    toc_total = time.perf_counter()

    pTime = f"Total processing time was  {toc_total - tic_total:0.4f} seconds"
    print(pTime)

if __name__ == "__main__":
    # argument: indir -> get all nc files in this directory
    # argument: map -> the offlinemap file already prepared
    # argument: outdir -> directory where remapped files should go
    # my_parser = arg.ArgumentParser()
    # my_parser.add_argument("--month", type=int)
    # my_parser.add_argument("--year", type=int)
    # args = my_parser.parse_args()
    
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--month", type=int, default=1)
    my_parser.add_argument("--year", type=int, default=2010)
    my_parser.add_argument("--day", type=int, default=1)
    my_parser.add_argument("--hour", type=int, default=99)
    args = my_parser.parse_args()
    main(args.year, args.month, args.day, args.hour )
