#!/usr/bin/env python
# Import packages 
import sys
import argparse as arg
# import modules in other directories
sys.path.append('/glade/work/juliob/PyRegridding/Regridder/')
sys.path.append('/glade/work/juliob/PyRegridding/Utils/')


import importlib
import glob
import copy
#import time
import os 
import subprocess as sp

import xarray as xr
import numpy as np
import pandas as pd

import ESMF as E

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



def Hregrid(case,BaseDir,ymdPat='*'):
    

    hsPat = 'cam.h0'
    
    ####################
    #case = "test2.04"
    #BaseDir = "/glade/derecho/scratch/juliob/archive/"


    #####################
    Src     = 'ne30pg3'
    SrcDir  = BaseDir+case+'/atm/hist/'
    SrcFile = SrcDir + case + '.cam.h0.2000-01.nc'

    ####################
    Dst     = 'fv0.9x1.25'
    DstTag  = '/regridded/'   # new subdirectory tag for regridded data
    DstDir = SrcDir.replace('/hist/', DstTag )

    #DstDir  = BaseDir+case+'/atm/regridded/'
    #DstFile = DstDir + case + '.cam.h0.2000-01.nc'

    
    #######
    os.makedirs( DstDir , exist_ok=True )
    
    ####################
    if (Src == 'ne30pg3'):
        srcHkey = 'c'
        src_type='mesh'
        src_scrip = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/ne30pg3_scrip_170611.nc'
        src_TopoFile = '/glade/p/cgd/amp/juliob/bndtopo/latest/ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_20230105.nc'

    if ((Dst == 'fv0.9x1.25') or (Dst=='fv1x1')):
        dstHkey = 'yx'
        dst_type='grid'
        dst_scrip = '/glade/p/cesmdata/cseg/inputdata/share/scripgrids/fv0.9x1.25_141008.nc'
        dst_TopoFile = '/glade/p/cesmdata/cseg/inputdata/atm/cam/topo/fv_0.9x1.25_nc3000_Nsw042_Nrs008_Co060_Fi001_ZR_160505.nc'

    # ----------------------------------------------
    # Make object for ESMF regridding from SRC
    # grid to CAM target. Scrip files need to be provided even 
    # when a weight file is used.
    # We make both conservative and bilinear mapping
    # files, because we will use conservative for 2D surface vars, but
    # bilinear for 3D met fields.
    # ----------------------------------------------
    RegridMethod = "BILINEAR"
    regrdB, srcfB, dstfB = erg.Regrid( srcScrip = src_scrip , 
                                    srcType  = src_type  ,
                                    dstScrip = dst_scrip ,
                                    dstType  = dst_type  ,
                                    RegridMethod = RegridMethod )

    RegridMethod = "CONSERVE"
    regrdC, srcfC, dstfC = erg.Regrid( srcScrip = src_scrip , 
                                    srcType  = src_type  ,
                                    dstScrip = dst_scrip ,
                                    dstType  = dst_type  ,
                                    RegridMethod = RegridMethod )

    ######################################
    # Klugy way to get lats and lons for
    # destination FV grid. SHould use
    # scrip files somehow
    ######################################
    DstTopoData  = xr.open_dataset( dst_TopoFile )
    
    
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


    ######################################
    # Get invariant grid stuff for FV
    ######################################
    lon_Dst = DstTopoData['lon'].values
    lat_Dst = DstTopoData['lat'].values

    slon=(lon_Dst[1:]+lon_Dst[0:-1] )/2.
    slat=(lat_Dst[1:]+lat_Dst[0:-1] )/2.

    lev = SrcData['lev'].values
    ilev = SrcData['ilev'].values

    nt,nz,ny,nx = 1,len(lev),len(lat_Dst),len(lon_Dst)
    print( nt,nz,ny,nx )

    dims   = ["lon","lat","time","lev","ilev","nbnd"]

    #########################################################
    # Lis of variables to be regridded.
    # This list should be a super-set of the one in the ADF
    # config..YAML.
    # PS needs to be here fro 3D vars.
    #########################################################
    ilist = [ 'SWCF'
            , 'PRECT'
            , 'LWCF'
            , 'PRECC'
            , 'PRECL'
            , 'PSL'
            , 'PS'
            , 'Q'
            , 'U'
            , 'T'
            , 'RELHUM'
            , 'TREFHT'
            , 'TS'
            , 'TAUX'
            , 'TAUY'
            , 'FSNT'
            , 'FLNT'
              , 'TMQ'
              , 'Nudge_U'
              , 'Nudge_V'
              , 'Nudge_T'
              , 'UTGW_MOVMTN'
              , 'VTGW_MOVMTN'
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

            DstData['time_bnds'] = SrcData['time_bnds']
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
                    # Conservative remapping for 2D vars
                    if (len_shap == 2 ):    
                        regrd = regrdC 
                        srcF  = srcfC 
                        dstF  = dstfC 
                        xfld_Dst = np.zeros( (nt,ny,nx) , dtype=np.float64 )
                        nlev=1
                        dims = ('time','lat','lon',)
                        Slice_Src = xfld_Src[0,:]
                        Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                            regrd = regrd , 
                                            srcField= srcF , 
                                            dstField= dstF , 
                                            srcGridkey= srcHkey ,
                                            dstGridkey= dstHkey )
                        xfld_Dst[0,:,:] = Slice_Dst

                    #############
                    # Bilinear remapping for 3D vars
                    if (len_shap == 3 ):    
                        regrd = regrdB 
                        srcF  = srcfB
                        dstF  = dstfB 
                        xfld_Dst = np.zeros( (nt,nz,ny,nx) , dtype=np.float64 )
                        nlev = nz
                        dims = ('time','lev','lat','lon',)
                        for L in np.arange( nlev ):
                            Slice_Src = xfld_Src[0,L,:]
                            Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                                        regrd = regrd , 
                                                        srcField= srcF , 
                                                        dstField= dstF , 
                                                        srcGridkey= srcHkey ,
                                                        dstGridkey= dstHkey )
                            xfld_Dst[0,L,:,:] = Slice_Dst

                    Dar = xr.DataArray( data=xfld_Dst , 
                                        dims=dims,
                                        attrs=attrs ,) 
                                
                    DstData[ fld ]= Dar
                    
                else:                     
                    print( f"No field {fld} in Src file" )
                    
            """
            DstData.to_netcdf( DstFile  , format='NETCDF3_CLASSIC'  ) #,encoding={'time': {'unlimited': True}} )
            ####################################################
            # Have to resort to 'ncks' $&%^!! to make time an
            # 'UNLIMITED' dimension as required by nrcat $&%^!!
            # This is also the SLOWEST part of this process.
            ####################################################
            sp.run(f" ncks --mk_rec_dmn time {DstFile} -o {DstDir}tmp.nc && \
                      mv {DstDir}tmp.nc {DstFile}", shell=True )
            print(f"Made time dim UNLIMITED using ncks " )
            """
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
    
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--case",     type=str )
    my_parser.add_argument("--BaseDir",  type=str, default="/glade/derecho/scratch/juliob/archive/")
    my_parser.add_argument("--ymdPat",  type=str, default="*")
    args = my_parser.parse_args()
    Hregrid( case=args.case, BaseDir=args.BaseDir, ymdPat=args.ymdPat )
