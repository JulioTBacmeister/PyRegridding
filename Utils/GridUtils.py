import sys
workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} ")

import numpy as np
import xarray as xr

from PyRegridding.Utils import MyConstants as C



pi = C.pi()

#####################
def area2d(lon,lat):
    #inputs 
    #   lon 1D lon vector
    #   lat 1D lat vector
    
    nx = np.size( lon )
    ny = np.size( lat )
    area = np.zeros( (ny, nx) )
    latr = (pi/180.) * lat
    
    for j in np.arange( ny ):
        area[j,:] = np.cos(latr[j])
        
    area = ( 4*pi / np.sum(area) ) * area
    
    return area

########################
def scrip_etc(grid=None):
    #########################
    # scrip files and other 
    # grid info
    #########################
    if (grid == 'ne30pg3'):
        Hkey='c'
        gtype='mesh'
        scrip = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/ne30pg3_scrip_170611.nc'
    elif ((grid == 'fv0.9x1.25') or (grid=='fv1x1')):
        Hkey='yx'
        gtype='grid'
        scrip = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.9x1.25_141008.nc'
    elif ( (grid == 'mpas120') or (grid == 'mpasa120') ) :
        Hkey='c'
        gtype='mesh'
        scrip = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/mpasa120_SCRIP_desc_211008.nc'
    else:
        scrip=''
        Hkey=''
        gtype=''

    return scrip,Hkey,gtype

##############################
def latlon(grid=None, scrip=None,Hkey=None ,get_area=False):

    if ( grid is not None ):
        Ginfo = gridInfo( grid )
        scrip=Ginfo['scrip']
        Hkey=Ginfo['Hkey']
    
    S = xr.open_dataset( scrip )
    
    if (Hkey == 'c'):
        lon = S.grid_center_lon.values
        lat = S.grid_center_lat.values
        if (get_area == True ):
            area = S.grid_area.values
            return lat,lon,area
        else:
            return lat,lon
    
    if (Hkey == 'yx'):
        lon = np.unique( S.grid_center_lon.values )
        lat = np.unique( S.grid_center_lat.values )
        return lat,lon
        
        
        
    return lat,lon

#############################
def gridInfo( grid=None , **kwargs ):

    cesm_inputdata_dir = '/glade/campaign/cesm/cesmdata/cseg/inputdata/'
    myGridFiles = '/glade/work/juliob/GridFiles/'

    if ('Vgrid' in kwargs):
        Vgrid = kwargs['Vgrid']
        if (Vgrid == 'L135' ):
            # Read in CAM L135 vertical grid
            # This is the "healed" WACCM L135 grid, top at 140km
            VgridFile = f'{myGridFiles}/Vertical/GRID_135L_CAM7_OrigAB_c20241011.nc'
        if (Vgrid == 'L120' ):
            # Read in CAM L120 vertical grid
            # This is a truncated version of the "healed" WACCM L135 grid, top at 85km
            VgridFile = f'{myGridFiles}/Vertical/GRID_120L_CAM7_OrigAB_Truncated_L135_c20241011.nc'
        if (Vgrid == 'L94_truncated_110' ):
            # Read in CAM L94 vertical grid
            # Truncated version of WACCM L110
            VgridFile = f'{myGridFiles}/Vertical/GRID_94L_CAM7_WACCM_Truncated_L110_c20241016.nc'
        if (Vgrid == 'L93' ):
            # Read in CAM L93 vertical grid
            # The file below though poorly named is the correct one to use for L93
            VgridFile = f'{myGridFiles}/Vertical/GRID_93L_CAM7_OrigAB_c20240514.nc'
        if (Vgrid == 'L58' ):
            # Read in CAM L58 vertical grid
            VgridFile = f'{myGridFiles}/Vertical/GRID_48_taperstart10km_lowtop_BL10_v3p1_beta1p75.nc'
        if (Vgrid == 'L56_86km' ):
            # Read in a truncated WACCM (v6 70-level) grid with a top near 86km
            VgridFile = f'{myGridFiles}/Vertical/GRID_56L_CAM7_TruncatedWACCM-Top86km_c20240705.nc'
        if (Vgrid == 'L42' ):
            # Not sure what this 32 level+ PBL refinements ??? Check!! (July 5 2024)
            VgridFile = f'{myGridFiles}/Vertical/GRID_42L_CAM7.nc'
        if (Vgrid == 'L32' ):
            # 32-level CAM6 grid
            VgridFile = f'{myGridFiles}/Vertical/GRID_32L_CAM6.nc'
    else:
        VgridFile = ''

    if ('VgridOnly' in kwargs):
        return VgridFile
    
    # ---------------------------------------------
    # This if block is to account
    # for the fact that sometimes (most times)
    # when the destination grid is neXXnp4 or
    # we are making IC files for an neXXpg3 run. 
    # The topo that should be used for this is 
    # 'PHIS_gll' on the neXXpg3 TopoFile ...
    #----------------------------------------------
    if ('IC_for_pg' in kwargs):
        IC_for_pg_ = kwargs['IC_for_pg']
    else:
        IC_for_pg_ = False

    # Remember to make all conditions after first 'elif'!!!!
    #---------------------------------------------------------
    if (grid == 'ne16pg3'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne16pg3_scrip_170429.nc'
        TopoFile = 'N/A'
        p_00 = 100_000.

    elif (grid == 'ne30pg3'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne30pg3_scrip_170611.nc'
        #TopoFile = cesm_inputdata_dir+'atm/cam/topo/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_20230105.nc'
        TopoFile = '/glade/campaign/cgd/amp/pel/topo/files/se/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_noleak_20240117.nc'
        p_00 = 100_000.

    elif (grid == 'ne120pg3'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne120pg3_scrip_170628.nc'
        TopoFile = myGridFiles+'Topo/ne120pg3_gmted2010_modis_bedmachine_nc3000_Laplace0025_noleak_20240326.nc'
        p_00 = 100_000.

    elif (grid == 'ne120np4'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne120np4_pentagons_100310.nc'
        if (IC_for_pg_ ==True):
            TopoFile = myGridFiles+'Topo/ne120pg3_gmted2010_modis_bedmachine_nc3000_Laplace0025_noleak_20240326.nc'
            print( f' Grabbed pg3 TopoFile even though grid is {grid}' )
        else:
            TopoFile = myGridFiles+'Topo/ne120np4_gmted2010_modis_bedmachine_nc3000_Laplace0025_noleak_20240326.nc'
        p_00 = 100_000.

    elif (grid == 'ne240pg3'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne240pg3_scrip_170628.nc'
        TopoFile = myGridFiles+'Topo/ne240pg3_gmted2010_modis_bedmachine_nc3000_Laplace0012_noleak_20240329.nc'
        p_00 = 100_000.

    elif (grid == 'ne240np4'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne240np4_091227_pentagons.nc'
        if (IC_for_pg_ ==True):
            TopoFile = myGridFiles+'Topo/ne240pg3_gmted2010_modis_bedmachine_nc3000_Laplace0012_noleak_20240329.nc'
            print( f' Grabbed pg3 TopoFile even though grid is {grid}' )
        else:
            TopoFile = 'N/A' #myGridFiles+'Topo/ne120np4_gmted2010_modis_bedmachine_nc3000_Laplace0025_noleak_20240326.nc'
        p_00 = 100_000.

    elif (grid == 'ne480pg3'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne480pg3_scrip_200108.nc'
        TopoFile = myGridFiles+'Topo/ne480pg3_gmted2010_modis_bedmachine_nc3000_NoAniso_Laplace0006_noleak_20240412.nc'
        p_00 = 100_000.

    elif (grid == 'ne480np4'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne480np4_scrip_c200409.nc'
        if (IC_for_pg_ ==True):
            TopoFile = myGridFiles+'Topo/ne480pg3_gmted2010_modis_bedmachine_nc3000_NoAniso_Laplace0006_noleak_20240412.nc'
            print( f' Grabbed pg3 TopoFile even though grid is {grid}' )
        else:
            TopoFile = 'N/A' #myGridFiles+'Topo/ne120np4_gmted2010_modis_bedmachine_nc3000_Laplace0025_noleak_20240326.nc'
        p_00 = 100_000.

    elif (grid == 'ne30np4'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne30np4_091226_pentagons.nc'
        if (IC_for_pg_ ==True):
            TopoFile = '/glade/campaign/cgd/amp/pel/topo/files/se/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_noleak_20240117.nc'
            print( f' Grabbed pg3 TopoFile even though grid is {grid}' )
        else:
            TopoFile = 'N/A' 
        p_00 = 100_000.

    elif (grid == 'ne4np4'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne4np4_pentagons.nc'
        TopoFile = 'N/A' #
        p_00 = 100_000.

    elif (grid == 'POLARRES'):
        Hkey = 'c'
        type ='mesh'
        scrip ='/glade/work/aherring/grids/var-res/ne0np4.POLARRES.ne30x4/grids/POLARRES_ne30x4_np4_SCRIP.nc'
        TopoFile = '/glade/work/aherring/grids/var-res/ne0np4.POLARRES.ne30x4/topo/POLARRES_gmted2010_modis_bedmachine_nc3000_Laplace0100_noleak_20240118.nc'
        p_00 = 100_000.
        
    elif (grid == 'Arctic'):
        Hkey = 'c'
        type='mesh'
        scrip = '/glade/work/aherring/grids/var-res/ne0np4.ARCTIC.ne30x4/grids/ne0ARCTICne30x4_scrip_191212.nc'
        TopoFile = cesm_inputdata_dir+'atm/cam/topo/se/ne30x4_ARCTIC_nc3000_Co060_Fi001_MulG_PF_RR_Nsw042_c200428.nc'
        p_00 = 100_000.

    elif (grid == 'MESO01'):
        Hkey = 'c'
        type='mesh'
        scrip = myGridFiles+'/Scrip/MESO01_ne30x4_np4_SCRIP.nc'
        TopoFile = myGridFiles+'/Topo/topo_ne0np4.MESO01.ne30x4_gmted2010_modis_bedmachine_nc3000_Laplace0100_noleak_230725.nc'
        p_00 = 100_000.

    elif (grid == 'MESO03'):
        Hkey = 'c'
        type='mesh'
        scrip = myGridFiles+'/Scrip/MESO03_ne30x4_np4_SCRIP.nc'
        TopoFile = 'N/A'
        p_00 = 100_000.

    elif ((grid == 'fv0.9x1.25') or (grid=='fv1x1')):
        Hkey = 'yx'
        type='grid'
        scrip = cesm_inputdata_dir+'share/scripgrids/fv0.9x1.25_141008.nc'
        TopoFile = cesm_inputdata_dir+'atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_160505.nc'
        p_00 = 100_000.

    elif ((grid == 'fv0.23x0.31') or (grid == 'fvQxQ') or (grid=='quarter degree') or (grid=='25km') ):
        # This grid exists primarily for plotting. So actual topo is not really needed
        Hkey = 'yx'
        type='grid'
        scrip = cesm_inputdata_dir+'share/scripgrids/fv0.23x0.31_071004.nc'
        TopoFile = 'N/A' #cesm_inputdata_dir+'atm/cam/topo/  ??? '
        p_00 = 100_000.

    elif ((grid == 'latlonOxO') or (grid=='eighth degree') or (grid=='12km') or (grid=='14km') ):
        # This grid exists primarily for plotting. So actual topo is not really needed
        Hkey = 'yx'
        type='grid'
        scrip = '/glade/work/juliob/GridFiles/Scrip/latlon_OxO_scrip.nc'
        TopoFile = 'N/A' #cesm_inputdata_dir+'atm/cam/topo/  ??? '
        p_00 = 100_000.

    elif (grid == 'CCIASI'  ):
        # This grid exists for CCIASI analysis. So actual topo is not really needed
        Hkey = 'yx'
        type='grid'
        scrip = '/glade/work/juliob/GridFiles/Scrip/CCIASI_scrip.nc'
        TopoFile = 'N/A' #cesm_inputdata_dir+'atm/cam/topo/  ??? '
        p_00 = 100_000.

    elif (grid == 'ERA5'):
        Hkey = 'yx'
        type='grid'
        scrip = '/glade/work/juliob/ERA5-proc/ERA5interp/grids/ERA5_640x1280_scrip.nc'
        TopoFile = '/glade/work/juliob/ERA5-proc/ERA5interp/phis/ERA5_phis.nc'
        p_00 = 1.0

    elif (grid == 'ERAI'):
        Hkey = 'yx'
        type='grid'
        scrip = '/glade/work/juliob/ERA-I-grids/ERAI_256x512_scrip.nc'
        TopoFile = '/glade/scratch/juliob/erai_2017/ei.oper.an.ml.regn128sc.2017010100.nc'
        p_00 = 100_000.
        
    else:
        Hkey = ''
        type=''
        scrip = ''
        TopoFile = ''
        p_00 = 100_000.

    
    Res = { 'Hkey' : Hkey ,
            'type' :  type ,
            'scrip' : scrip ,
            'TopoFile' : TopoFile ,
            'VgridFile' : VgridFile, 
            'p_00' : p_00 }
             
    return Res 
    
#############################
def gridKey( Var ):
    
    print( " IN working version ")
    VarDims = Var.dims
    ndim = len(VarDims)
    
    if (VarDims[0]=='time'):
        gridKey = 't'
    elif (VarDims[0]=='lev'):
        gridKey = 'z'
    elif (VarDims[0]=='lat'):
        gridKey = 'y'
    elif (VarDims[0]=='ncol'):
        gridKey = 'c'
    else:
        gridKey = 'not found'
        return gridKey
    
    if (ndim>1):
        if (VarDims[1]=='time'):
            gridKey = gridKey + 't'
        elif (VarDims[1]=='lev'):
            gridKey = gridKey + 'z'
        elif (VarDims[1]=='lat'):
            gridKey = gridKey + 'y'
        elif (VarDims[1]=='ncol'):
            gridKey = gridKey + 'c'
        else:
            gridKey = gridKey + ' not found'
            return gridKey
    
    if (ndim>2):
        if (VarDims[2]=='time'):
            gridKey = gridKey + 't'
        elif (VarDims[2]=='lev'):
            gridKey = gridKey + 'z'
        elif (VarDims[2]=='lat'):
            gridKey = gridKey + 'y'
        elif (VarDims[2]=='ncol'):
            gridKey = gridKey + 'c'
        else:
            gridKey = gridKey + ' not found'
            return gridKey
    
    if (ndim>3):
        if (VarDims[3]=='time'):
            gridKey = gridKey + 't'
        elif (VarDims[3]=='lev'):
            gridKey = gridKey + 'z'
        elif (VarDims[3]=='lat'):
            gridKey = gridKey + 'y'
        elif (VarDims[3]=='lon'):
            gridKey = gridKey + 'x'
        elif (VarDims[3]=='ncol'):
            gridKey = gridKey + 'c'
        else:
            gridKey = gridKey + ' not found'
            return gridKey
    
    
    return gridKey
