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

try:
    import ESMF as E
except ImportError:
    import esmpy as E

import scripGen as SG
import esmfRegrid as erg


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



def Hregrid(case,BaseDir,Dst,Src,ymdPat,hsPat,DstSubDir,clean=False, AllConservative=False ):
    

    #######################
    SrcDir  = BaseDir+case+'/atm/hist/'
    SrcFile = SrcDir + case + '.cam.h0.2000-01.nc'

    """
    ####################
    if ((Dst == 'fv0.9x1.25') or (Dst =='fv1x1')):
        DstTag  = '/regridded/'   # new subdirectory tag for regridded data
    ####################
    elif (Dst == 'fv0.23x0.31'):
        DstTag  = '/regridded_0.25/'   # new subdirectory tag for regridded data
    ####################
    else: 
        DstTag  = '/regridded/'   # new subdirectory tag for regridded data
    """
    DstTag = f'/{DstSubDir}/'
    DstDir = SrcDir.replace('/hist/', DstTag )

    #######
    os.makedirs( DstDir , exist_ok=True )


    if (clean==True ):
        # List all files in the directory
        for filename in os.listdir(DstDir):
            file_path = os.path.join(DstDir, filename)
            # Check if it's a file
            if os.path.isfile(file_path):
                # Delete the file
                os.remove(file_path)
            
    ############################################################
    # dst_TopoFile 's here are only used to get lats and lons
    # A better method may exist - GrU.latlon
    ############################################################
    DstInfo = GrU.gridInfo(Dst) #,Vgrid=DstVgrid)
    dstHkey = DstInfo['Hkey']
    dst_type =DstInfo['type']
    dst_scrip =DstInfo['scrip']
    dst_TopoFile = DstInfo['TopoFile']

    SrcInfo = GrU.gridInfo(Src)
    srcHkey = SrcInfo['Hkey']
    src_type =SrcInfo['type']
    src_scrip =SrcInfo['scrip']
    src_TopoFile = SrcInfo['TopoFile']
    print( f"Used NEW, concise gridInfo function .... ...." )

    # ----------------------------------------------
    # Make object for ESMF regridding from SRC
    # grid to CAM target. Scrip files need to be provided even 
    # when a weight file is used.
    # We make both conservative and bilinear mapping
    # files, because we will use conservative for 2D surface vars, but
    # bilinear for 3D met fields.
    # ----------------------------------------------
    print(f'Creating Bilinear map {src_scrip} => {dst_scrip}')
    RegridMethod = "BILINEAR"
    regrdB, srcfB, dstfB = erg.Regrid( srcScrip = src_scrip , 
                                    srcType  = src_type  ,
                                    dstScrip = dst_scrip ,
                                    dstType  = dst_type  ,
                                    RegridMethod = RegridMethod )

    print(f'Creating Conservative map {src_scrip} => {dst_scrip}')
    RegridMethod = "CONSERVE"
    regrdC, srcfC, dstfC = erg.Regrid( srcScrip = src_scrip , 
                                    srcType  = src_type  ,
                                    dstScrip = dst_scrip ,
                                    dstType  = dst_type  ,
                                    RegridMethod = RegridMethod )

    
    ###########################################################
    # Pattern to match all h0 files in the specified directory
    ###########################################################
    #SrcFile = SrcDir + case + '.cam.h0.*.nc'

    
    #SrcFile = SrcDir + case + '.cam.h0.' + ymdPat + '.nc'
    #DstFile = DstDir + case + '.cam.h0.' + ymdPat + '.nc'
    SrcFile = SrcDir + case + '.'+ hsPat + '.' + ymdPat + '.nc'
    DstFile = DstDir + case + '.'+ hsPat + '.' + ymdPat + '.nc'

    print(f" File list made from pattern: {SrcFile} " )
    
    # Get a list of all matching file paths
    file_list = sorted( glob.glob(SrcFile)  )
    dst_file_list = sorted( glob.glob(DstFile)  )

    ######################################
    # Open first dataset to get lev, ilev
    ######################################
    SrcData = xr.open_dataset( file_list[0] )
    #--------------------------------------
    # Is time a dimension. If not add a stub dimension
    # of length 1.
    #---------------------------------
    if ( 'time' not in SrcData.dims ):
        nt_Src=1
        pad_w_time=True
        # Add a time dimension with length 1
        SrcData = SrcData.expand_dims(time=1)
        # Create a time coordinate as a NumPy array
        time_coord = np.array([0.])

        # Assign the time coordinate to the dataset
        SrcData = SrcData.assign_coords(time=time_coord)
    else:
        nt_Src = SrcData.dims['time']
        pad_w_time=False
    
    ######################################
    # Get invariant grid stuff for FV
    ######################################
    if (dst_TopoFile != 'N/A' ):
        ######################################
        # Klugy way to get lats and lons for
        # destination FV grid. SHould use
        # scrip files somehow
        ######################################
        DstTopoData  = xr.open_dataset( dst_TopoFile )    
        lon_Dst = DstTopoData['lon'].values
        lat_Dst = DstTopoData['lat'].values
    else:
        lat_Dst,lon_Dst = GrU.latlon( scrip=dst_scrip, Hkey=dstHkey )
        print( f"Used GrU.latlon {dst_scrip} for lat lon ")

    slon=(lon_Dst[1:]+lon_Dst[0:-1] )/2.
    slat=(lat_Dst[1:]+lat_Dst[0:-1] )/2.

    lev = SrcData['lev'].values
    ilev = SrcData['ilev'].values

    nt,nz,ny,nx = nt_Src,len(lev),len(lat_Dst),len(lon_Dst)
    print( nt,nz,ny,nx )

    dims   = ["lon","lat","time","lev","ilev","nbnd"]

    #########################################################
    # List of variables to be regridded.
    # This list should be a super-set of the one in the ADF
    # config..YAML.
    # PS needs to be here for 3D vars.
    #########################################################
    ilist = [ 'area'
            , 'SWCF'
            , 'PRECT'
            , 'LWCF'
            , 'PRECC'
            , 'PRECL'
            , 'PSL'
            , 'PS'
            , 'Q'
            , 'U'
            , 'V'
            , 'T'
            , 'OMEGA'
            , 'OMEGAU'
            , 'OMEGAV'
            , 'RELHUM'
            , 'TREFHT'
            , 'TS'
            , 'TAUX'
            , 'TAUY'
            , 'FSNT'
            , 'FLNT'
              ,'SHFLX'
              ,'LHFLX'
              , 'TMQ'
              , 'Nudge_U'
              , 'Nudge_V'
              , 'Nudge_T'
              , 'UTEND_CORE'
              , 'VTEND_CORE'
              , 'UTEND_PHYSTOT'
              , 'VTEND_PHYSTOT'
              , 'UTEND_GWDTOT'
              , 'VTEND_GWDTOT'
              , 'UTGW_MOVMTN'
              , 'VTGW_MOVMTN'
              , 'UTGWORO'
              , 'UTGWSPEC'
              , 'BUTGWSPEC'
              , 'UPWP_CLUBB'
              , 'VPWP_CLUBB'
              , 'WP2_CLUBB'
              , 'WP3_CLUBB'
              , 'WPTHLP_CLUBB'
              , 'THLP2_CLUBB'
              , 'STEND_CLUBB'
              , 'PHIS'
              , 'U10'
              , 'LANDFRAC' ]


    print( file_list )

    ############################################
    # March through files to be horz regridded
    ############################################
    for fileN in file_list:
        # Set source file
        SrcFile = fileN
        
        # Make destination file name by replacing 'hist' in path
        DstFile = SrcFile.replace('/hist/', DstTag )

        DstFileSizeThresh = 0
        """
        if (DstFile in dst_file_list):
            SkipThisFile = True
        else:
            SkipThisFile = False
        """
        DstFileAlreadyExists = os.path.exists(DstFile)
        if (DstFileAlreadyExists==True):
            DstFileSize = 1000 # os.path.getsize(DstFile)
        else:
            DstFileSize = 0
        
        if ( (DstFileAlreadyExists==True) and (DstFileSize > DstFileSizeThresh) ):
            SkipThisFile = True
        else:
            SkipThisFile = False
            
        
        if (SkipThisFile == True ):
            print( f" {DstFile} Already exists " )
        else:
            SrcData = xr.open_dataset( SrcFile )
            print( f"Read {SrcFile}")

            time = SrcData['time'].values
            timedim = SrcData['time'].dims

            coords = dict( 
                time = ( ["time"],  time ),
                lon  = ( ["lon"],lon_Dst),
                lat  = ( ["lat"],lat_Dst ),
                slon  = ( ["slon"],slon ),
                slat  = ( ["slat"],slat ),
                lev  = ( ["lev"],lev),
                ilev = ( ["ilev"],ilev),
                nbnd = ( ["nbnd"], np.array([0,1] ) ),
            )


            DstData = xr.Dataset( coords=coords  )

            if ( 'time_bnds' in SrcData ):
                DstData['time_bnds'] = SrcData['time_bnds']
            if ( 'time_bounds' in SrcData ):
                DstData['time_bounds'] = SrcData['time_bounds']
            DstData['date'] = SrcData['date']
            DstData['datesec'] = SrcData['datesec']

            DstData['hyai'] = SrcData['hyai']
            DstData['hybi'] = SrcData['hybi']
            DstData['hyam'] = SrcData['hyam']
            DstData['hybm'] = SrcData['hybm']

            if ( ('PRECT' not in SrcData) and ( ('PRECL' in SrcData) and ('PRECC' in SrcData)   ) ):
                SrcData['PRECT'] = SrcData['PRECC'] + SrcData['PRECL']
                new_attributes = {
                    'units': 'ms-1',
                    'description': 'sum of PRECC and PRECL',}
                SrcData['PRECT'].attrs.update(new_attributes)
                print(f"Created PRECT")
            if ( ('OMEGAU' not in SrcData) and ( ('OMEGA' in SrcData) and ('U' in SrcData)   ) ):
                SrcData['OMEGAU'] = SrcData['OMEGA'] * SrcData['U']
                new_attributes = {
                    'units': 'Pa ms-2',
                    'description': 'Vertical flux of zonal mom.',}
                SrcData['OMEGAU'].attrs.update(new_attributes)
                print(f"Created OMEGAU")

            ###################
            # March through variables
            ###################
            for fld in ilist:
                if (fld in SrcData):
                    print( f"Regridding {fld}" )
                    
                    shap = np.shape( SrcData[fld] ) 
                    len_shap = len( shap )
                    xfld_Src = SrcData[fld].values
                    attrs =  SrcData[fld].attrs 

                    ####################
                    # The following 2D vs 3D determination
                    # should be done in a better way.
                    # What is here will work for SE CAM output, i.e., 
                    # (time,col) or (time,lev,col)
                    #######################
                    # Always Conservative remapping for 2D vars and area
                    if (fld == 'area' ):    
                        regrd = regrdC 
                        srcF  = srcfC 
                        dstF  = dstfC 
                        xfld_Dst = np.zeros( (ny,nx) , dtype=np.float64 )
                        nlev=1
                        dims = ('lat','lon',)
                        Slice_Src = xfld_Src[:]
                        Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                            regrd = regrd , 
                                            srcField= srcF , 
                                            dstField= dstF , 
                                            srcGridkey= srcHkey ,
                                            dstGridkey= dstHkey )
                        xfld_Dst[:,:] = Slice_Dst
                    if (len_shap == 2 ):    
                        regrd = regrdC 
                        srcF  = srcfC 
                        dstF  = dstfC 
                        xfld_Dst = np.zeros( (nt,ny,nx) , dtype=np.float64 )
                        nlev=1
                        dims = ('time','lat','lon',)
                        for tin in np.arange( nt ):
                            Slice_Src = xfld_Src[tin,:]
                            Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                                regrd = regrd , 
                                                srcField= srcF , 
                                                dstField= dstF , 
                                                srcGridkey= srcHkey ,
                                                dstGridkey= dstHkey )
                            xfld_Dst[tin,:,:] = Slice_Dst


                    #############
                    # Bilinear remapping for most 3D vars if AllConservative==False
                    if (len_shap == 3 ): 
                        if ( (fld not in ('OMEGAU','OMEGAV' )) and (AllConservative==False) ):
                            regrd = regrdB 
                            srcF  = srcfB
                            dstF  = dstfB 
                        else: 
                            regrd = regrdC 
                            srcF  = srcfC
                            dstF  = dstfC 
                            print(f" using conervative remapping for {fld}")
                            
                        xfld_Dst = np.zeros( (nt,nz,ny,nx) , dtype=np.float64 )
                        nlev = nz
                        dims = ('time','lev','lat','lon',)
                        for tin in np.arange( nt ):
                            for L in np.arange( nlev ):
                                Slice_Src = xfld_Src[tin,L,:]
                                Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                                            regrd = regrd , 
                                                            srcField= srcF , 
                                                            dstField= dstF , 
                                                            srcGridkey= srcHkey ,
                                                            dstGridkey= dstHkey )
                                xfld_Dst[tin,L,:,:] = Slice_Dst

                    Dar = xr.DataArray( data=xfld_Dst , 
                                        dims=dims,
                                        attrs=attrs ,) 
                                
                    DstData[ fld ]= Dar
                    
                else:                     
                    print( f"No field {fld} in Src file" )
                    
            #######################################################
            # ChatGPT says you can do this ... hallucinating???
            # When saving, specify that 'time' should be an unlimited dimension
            # ds_expanded.to_netcdf('outfile.nc', unlimited_dims='time')
            #
            # ChatGPT seems to be right. This goes through and 'ncdump -h DstFile' has
            #
            # dimensions:
            #   time = UNLIMITED ; // (1 currently)
            #
            # We'll see if ADF takes this
            # Does just fine .... YAY! THANKS CHATGPT !!!! You're
            # the greatest!
            ########################################################
            DstData.to_netcdf( DstFile  , unlimited_dims='time' )
            print(f"Made time dim UNLIMITED using xarray to_netcdf  " )
            

        print( f"Done {DstFile}" )

if __name__ == "__main__":
    
    #####################
    # Src     = 'ne30pg3'
    # Dst     = 'fv0.9x1.25'
    # hsPat = 'cam.h0'
    


    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--case",     type=str )
    my_parser.add_argument("--BaseDir",  type=str, default="/glade/derecho/scratch/juliob/archive/")
    my_parser.add_argument("--ymdPat",  type=str, default="*")
    my_parser.add_argument("--hsPat",  type=str, default="cam.h0")
    my_parser.add_argument("--Src",  type=str, default="ne30pg3")
    my_parser.add_argument("--Dst",  type=str, default="fv0.9x1.25")
    my_parser.add_argument("--DstSubDir", type=str, default="regridded")
    my_parser.add_argument("--clean", action="store_true", help="Set clean to True when --clean")
    my_parser.add_argument("--AllConservative", action="store_true", help="Set AllConservative to True when --AllConservative")

    args = my_parser.parse_args()
    Hregrid( case=args.case, BaseDir=args.BaseDir, ymdPat=args.ymdPat, hsPat=args.hsPat, Src=args.Src, Dst=args.Dst, DstSubDir=args.DstSubDir, 
            clean=args.clean, AllConservative=args.AllConservative )
