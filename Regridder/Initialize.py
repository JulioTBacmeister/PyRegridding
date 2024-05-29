#!/usr/bin/env python
# Import packages 
import sys

# This contains all the data
#----------------------------
import GlobalVarClass
from GlobalVarClass import Gv

import xarray as xr
import numpy as np

try:
    import ESMF as E
except ImportError:
    import esmpy as E

import importlib
import time

import esmfRegrid as erg


sys.path.append('../Utils/')
import GridUtils as GrU

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
#
# Note that in this code you should regard 'ERA' as the 'source'
# and 'CAM' as the 'destination'.  This is inherited and should
# be cleaned up
#-------------------------------------------------------------

def prep(Dst = 'ne30pg3', DstVgrid='L58',  Src='ERA5', WOsrf=False , RegridMethod="CONSERVE", IC_for_pg=False ):
    #---------------------------------------------
    # This function sets-up variables and objects 
    # that are need for horizontal and vertical 
    # regridding of ERA reanalyses.
    #---------------------------------------------    
    #------- 
    # Begin
    #-------
    
    
    tic_overall = time.perf_counter()
    Gv.MyDst,Gv.MyDstVgrid,Gv.MySrc = Dst,DstVgrid,Src

    Gv.RegridMethod = RegridMethod
    Gv.doWilliamsonOlson = WOsrf
    Gv.p_00_CAM = 100_000.

    cesm_inputdata_dir = '/glade/campaign/cesm/cesmdata/cseg/inputdata/'
    #my_bndtopo = 
    print( f"In prep Src= {Src} to Dst={Dst} " )

    DstInfo = GrU.gridInfo(Dst,Vgrid=DstVgrid, IC_for_pg= IC_for_pg )
    Gv.dstHkey = DstInfo['Hkey']
    Gv.dst_type =DstInfo['type']
    Gv.dst_scrip =DstInfo['scrip']
    Gv.dst_TopoFile = DstInfo['TopoFile']
    Gv.dstVgridFile = DstInfo['VgridFile' ] 

    SrcInfo = GrU.gridInfo(Src)
    Gv.srcHkey = SrcInfo['Hkey']
    Gv.src_type =SrcInfo['type']
    Gv.src_scrip =SrcInfo['scrip']
    Gv.src_TopoFile = SrcInfo['TopoFile']
    Gv.p_00_ERA = SrcInfo['p_00']
    print( f"Used NEW, concise gridInfo function .... ...." )

    # Set grid keys for Src ERA5 reanalysis
    Gv.srcTHkey  = 't'  + Gv.srcHkey
    Gv.srcZHkey  = 'z'  + Gv.srcHkey
    Gv.srcTZHkey = 'tz' + Gv.srcHkey

    # Set grid keys for Dst CAM-SE
    Gv.dstTHkey  = 't'  + Gv.dstHkey
    Gv.dstZHkey  = 'z'  + Gv.dstHkey
    Gv.dstTZHkey = 'tz' + Gv.dstHkey
 
    # ----------------------------------------------
    # Get all topo data we will use
    # Read in CAM topography. Also get
    # lon and lat and area for CAM (Dst)
    # grid.
    # 2024-03-26:
    # Here we need to account for the fact 
    # that neXXpg3 grids need ICs on GLL (neXXnp4) 
    # points. The right topo for generating these ICs
    # lives in the neXXpg3 topo file under a different name:
    #      PHIS_gll
    # ----------------------------------------------
    dsTopo_CAM=xr.open_dataset( Gv.dst_TopoFile )
    varsCAM  = list( dsTopo_CAM.variables )
    if (IC_for_pg==False):
        Gv.phis_CAM = dsTopo_CAM['PHIS'].values
    else:
        Gv.phis_CAM = dsTopo_CAM['PHIS_gll'].values
        print(f'   -- Making Initial Condition file for pg3 analog of {Dst}. Using PHIS_gll in TopoFile' )
    #---------------------------------------
    # It would be cleaner to get lat,lon directly
    # from the SCRIP file
    # 2024-03-26:
    # IMplementing the above for unstructured grids.
    # Not sure I trust np.unique for lat-lon grids. 
    # Need to figure out about area ....
    #---------------------------------------
    if ( Gv.dstHkey == 'c' ):
        Gv.lat_CAM, Gv.lon_CAM, Gv.area_CAM = GrU.latlon( scrip = Gv.dst_scrip , Hkey=Gv.dstHkey, get_area=True )
    else:
        Gv.lon_CAM  = dsTopo_CAM['lon'].values
        Gv.lat_CAM  = dsTopo_CAM['lat'].values
        if ('area' in varsCAM):
            Gv.area_CAM = dsTopo_CAM['area'].values
        else:
            Gv.area_CAM = GrU.area2d( lon=Gv.lon_CAM, lat=Gv.lat_CAM )

    if (Src == 'ERA5'):
        # Read in ERA5 topography
        dsTopo_ERA=xr.open_dataset( Gv.src_TopoFile )
        Gv.phis_ERA=dsTopo_ERA['Z_GDS4_SFC'].values

    if (Src == 'ERAI'):
        # Read in ERA-I topography
        dsTopo_ERA=xr.open_dataset( Gv.src_TopoFile )
        Gv.phis_ERA=dsTopo_ERA['Z_GDS4_HYBL'].values

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



    # ----------------------------------------------
    #  Set-up regridding machinery
    # ----------------------------------------------
    # Scrip file for ERA5 created by ERA5scrip.ipynb
    if (Src == 'ERA5'):
        dsERAscrip = xr.open_dataset( Gv.src_scrip )
        Gv.area_ERA = -9999. #np.reshape( dsERAscrip['grid_area'].values , np.shape( phis_ERA ) )
    else:
        Gv.area_ERA = -9999.
        
    # ----------------------------------------------
    # Make object for ESMF regridding from SRC
    # grid to CAM target. Scrip files need to be provided even 
    # when a weight file is used
    # ----------------------------------------------
    Gv.regrd, Gv.srcf, Gv.dstf = erg.Regrid( srcScrip = Gv.src_scrip , 
                                    srcType  = Gv.src_type  ,
                                    dstScrip = Gv.dst_scrip ,
                                    dstType  = Gv.dst_type  ,
                                    write_weights = write_weights ,
                                    read_weights = read_weights ,
                                    weights_file = wgts_file_Con ,
                                    RegridMethod = RegridMethod )
    


    vCAM=xr.open_dataset( Gv.dstVgridFile )
    Gv.amid_CAM = vCAM['hyam'].values
    Gv.bmid_CAM = vCAM['hybm'].values
    Gv.aint_CAM = vCAM['hyai'].values
    Gv.bint_CAM = vCAM['hybi'].values

    print( f" Src scripfile {Gv.src_scrip} " )
    print( f" Dst scripfile {Gv.dst_scrip} " )
    print( f" Src topo file {Gv.src_TopoFile} " )
    print( f" Dst topo file {Gv.dst_TopoFile} " )
    print( f" {DstVgrid} Dst vertical grid from {Gv.dstVgridFile} " )


    toc = time.perf_counter()
    pTime = f"Prepping for {Src} to {Dst} proc in {__name__} took  {toc - tic_overall:0.4f} seconds"
    print(pTime)
 
    code = 1
    return code

