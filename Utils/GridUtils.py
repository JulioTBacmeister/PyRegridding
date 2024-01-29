import numpy as np
import MyConstants as C
import xarray as xr

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
def latlon(scrip,Hkey):
    
    S = xr.open_dataset( scrip )
    
    if (Hkey == 'c'):
        lon = S.grid_center_lon.values
        lat = S.grid_center_lat.values
        
    if (Hkey == 'yx'):
        lon = np.unique( S.grid_center_lon.values )
        lat = np.unique( S.grid_center_lat.values )
        
        
        
    return lat,lon

#############################
def gridInfo( grid=None , **kwargs ):

    cesm_inputdata_dir = '/glade/campaign/cesm/cesmdata/cseg/inputdata/'

    if (grid == 'ne30pg3'):
        Hkey = 'c'
        type='mesh'
        scrip = cesm_inputdata_dir+'share/scripgrids/ne30pg3_scrip_170611.nc'
        TopoFile = cesm_inputdata_dir+'atm/cam/topo/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_20230105.nc'
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

    elif ((grid == 'fv0.9x1.25') or (grid=='fv1x1')):
        Hkey = 'yx'
        type='grid'
        scrip = cesm_inputdata_dir+'share/scripgrids/fv0.9x1.25_141008.nc'
        TopoFile = cesm_inputdata_dir+'atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_160505.nc'
        p_00 = 100_000.

    elif (grid == 'fv0.23x0.31'):
        Hkey = 'yx'
        type='grid'
        scrip = cesm_inputdata_dir+'share/scripgrids/fv0.23x0.31_071004.nc'
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

    if ('Vgrid' in kwargs):
        Vgrid = kwargs['Vgrid']
        if (Vgrid == 'L93' ):
            # Read in CAM L93 vertical grid
            VgridFile = '/glade/work/juliob/ERA5-proc/CAM-grids/Vertical/GRID_93L_CAM7_c202312.nc'
        if (Vgrid == 'L58' ):
            # Read in CAM L58 vertical grid
            VgridFile = '/glade/work/juliob/ERA5-proc/CAM-grids/Vertical/GRID_48_taperstart10km_lowtop_BL10_v3p1_beta1p75.nc'
        if (Vgrid == 'L32' ):
            #Gv.dstVgridFile = cesm_inputdata_dir+'atm/cam/inic/se/f.e22.FC2010climo.ne30pg3_ne30pg3_mg17.cam6_2_022.002.cam.i.0020-01-01-00000_c200610.nc'
            VgridFile = '/glade/work/juliob/ERA5-proc/CAM-grids/Vertical/GRID_32L_CAM6.nc'
    else:
        VgridFile = ''
    
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