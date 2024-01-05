#!/usr/bin/env python
# Import packages 
import sys

# This contains all the data
#----------------------------
import GlobalVarClass
from GlobalVarClass import Gv

import xarray as xr
import numpy as np

import ESMF as E

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

def prep(Dst = 'ne30pg3', DstVgrid='L58',  Src='ERA5', WOsrf=False , RegridMethod="CONSERVE" ):
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

    Gv.doWilliamsonOlson = WOsrf
    Gv.p_00_CAM = 100_000.

    cesm_inputdata_dir = '/glade/campaign/cesm/cesmdata/cseg/inputdata/'
    #my_bndtopo = 
    print( f"In prep Src= {Src} to Dst={Dst} " )
    if (Dst == 'Arctic'):
        Gv.dstHkey = 'c'
        Gv.dst_type='mesh'
        Gv.dst_scrip = '/glade/work/aherring/grids/var-res/ne0np4.ARCTIC.ne30x4/grids/ne0ARCTICne30x4_scrip_191212.nc'
        Gv.dst_TopoFile = cesm_inputdata_dir+'atm/cam/topo/se/ne30x4_ARCTIC_nc3000_Co060_Fi001_MulG_PF_RR_Nsw042_c200428.nc'

    if (Dst == 'ne30pg3'):
        Gv.dstHkey = 'c'
        Gv.dst_type='mesh'
        Gv.dst_scrip = cesm_inputdata_dir+'share/scripgrids/ne30pg3_scrip_170611.nc'
        #Gv.dst_TopoFile = '/glade/p/cgd/amp/juliob/bndtopo/latest/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_20230105.nc'
        Gv.dst_TopoFile = cesm_inputdata_dir+'atm/cam/topo/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_20230105.nc'

    if ((Dst == 'fv0.9x1.25') or (Dst=='fv1x1')):
        Gv.dstHkey = 'yx'
        Gv.dst_type='grid'
        Gv.dst_scrip = cesm_inputdata_dir+'share/scripgrids/fv0.9x1.25_141008.nc'
        #dst_TopoFile='/glade/p/cgd/amp/juliob/bndtopo/latest/fv_0.9x1.25_gmted2010_modis_bedmachine_nc3000_Laplace0100_20220708.nc'
        Gv.dst_TopoFile = cesm_inputdata_dir+'atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_160505.nc'

    if (Src == 'ERA5'):
        Gv.srcHkey = 'yx'
        Gv.src_type='grid'
        Gv.src_scrip = '/glade/work/juliob/ERA5-proc/ERA5interp/grids/ERA5_640x1280_scrip.nc'
        Gv.src_TopoFile = '/glade/work/juliob/ERA5-proc/ERA5interp/phis/ERA5_phis.nc'
        Gv.p_00_ERA = 1.0

    if (Src == 'ERAI'):
        Gv.srcHkey = 'yx'
        Gv.src_type='grid'
        Gv.src_scrip = '/glade/work/juliob/ERA-I-grids/ERAI_256x512_scrip.nc'
        Gv.src_TopoFile = '/glade/scratch/juliob/erai_2017/ei.oper.an.ml.regn128sc.2017010100.nc'
        Gv.p_00_ERA = 100_000.

    # ----------------------------------------------
    # Get DST vertical grid from a file.
    # These should be small files but 
    # they aren't always.
    # ----------------------------------------------
    if (DstVgrid == 'L93' ):
        # Read in CAM L93 vertical grid
        Gv.dstVgridFile = '/glade/work/juliob/ERA5-proc/CAM-grids/Vertical/GRID_93L_CAM7_c202312.nc'

    if (DstVgrid == 'L58' ):
        # Read in CAM L58 vertical grid
        Gv.dstVgridFile = '/glade/work/juliob/ERA5-proc/CAM-grids/Vertical/GRID_48_taperstart10km_lowtop_BL10_v3p1_beta1p75.nc'

    if (DstVgrid == 'L32' ):
        #Gv.dstVgridFile = cesm_inputdata_dir+'atm/cam/inic/se/f.e22.FC2010climo.ne30pg3_ne30pg3_mg17.cam6_2_022.002.cam.i.0020-01-01-00000_c200610.nc'
        Gv.dstVgridFile = '/glade/work/juliob/ERA5-proc/CAM-grids/Vertical/GRID_32L_CAM6.nc'

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
    # ----------------------------------------------
    dsTopo_CAM=xr.open_dataset( Gv.dst_TopoFile )
    varsCAM  = list( dsTopo_CAM.variables )
    Gv.phis_CAM = dsTopo_CAM['PHIS'].values
    #---------------------------------------
    # It would be cleaner to get lat,lon directly
    # from the SCRIP file
    #---------------------------------------
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

