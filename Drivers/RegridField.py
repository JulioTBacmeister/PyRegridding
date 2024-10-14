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



def Horz(Dst,Src, xfld_Src=None, RegridMethod=None, RegridObj_In=None, RegridObj_Out=False ):
    
    ##############################################################
    #  We will refer to '2D' and '3D' fields regardless of the 
    #  dimensions of the array representing them. '2D' will simply 
    #  mean fields with no vertical or time dimension.
    ############################################################
    DstInfo = GrU.gridInfo(Dst) #,Vgrid=DstVgrid)
    dstHkey = DstInfo['Hkey']
    dst_type =DstInfo['type']
    dst_scrip =DstInfo['scrip']

    SrcInfo = GrU.gridInfo(Src)
    srcHkey = SrcInfo['Hkey']
    src_type =SrcInfo['type']
    src_scrip =SrcInfo['scrip']

        
    
    print( f"Used NEW, concise gridInfo function .... ...." )

    # ----------------------------------------------
    # Make object for ESMF regridding from SRC
    # grid to CAM target. Scrip files need to be provided even 
    # when a weight file is used.
    # We make both conservative and bilinear mapping
    # files, because we will use conservative for 2D surface vars, but
    # bilinear for 3D met fields.
    # ----------------------------------------------

    if ( RegridObj_In == None ):
        if (RegridMethod==None):
            RegridMethod_ = 'CONSERVE_2ND' 
        else:
            RegridMethod_ = RegridMethod
    
        regrd, srcf, dstf = erg.Regrid( srcScrip = src_scrip , 
                                        srcType  = src_type  ,
                                        dstScrip = dst_scrip ,
                                        dstType  = dst_type  ,
                                        RegridMethod = RegridMethod_ )
        if ( (RegridObj_Out == True) or (xfld_Src is None) ):
            print( f" Not interpolating. Returning: regrd, srcf, dstf " )
            return regrd, srcf, dstf
    else:
        print( f" Getting (regrd, srcf, dstf) from argument " ) 
        regrd, srcf, dstf = RegridObj_In[0], RegridObj_In[1], RegridObj_In[2]   

    src_shape = np.shape( xfld_Src )
    len_src_shape = len( src_shape )
    nzot,nz = 0,0
    if ( (len_src_shape == 1 ) and (srcHkey=='c' ) ):
        nzot = 0
        srcShape='c'
    elif ( (len_src_shape == 2) and (srcHkey=='c' ) ):
        nzot = src_shape[0]
        srcShape = 'oc'
    elif ( (len_src_shape == 2) and (srcHkey=='yx' ) ):
        nzot = 0
        srcShape = 'yx'
    elif ( (len_src_shape == 3) and (srcHkey=='c' ) ):
        nzot = src_shape[0]
        nz = src_shape[1]
        srcShape = 'tzc'
    elif ( (len_src_shape == 3) and (srcHkey=='yx' ) ):
        nzot = src_shape[0]
        srcShape = 'oyx'
    elif ( (len_src_shape == 4) and (srcHkey=='yx' ) ):
        nzot = src_shape[0]
        nz = src_shape[1]
        srcShape = 'tzyx'
    else:
        print("your array can't be classified " )
        srcShape='nowayjose'
    
    
    ###########################################################

    lat_Dst,lon_Dst = GrU.latlon( scrip=dst_scrip, Hkey=dstHkey )
    print( f"Used GrU.latlon {dst_scrip} for lat lon ")

    ncol,ny,nx=0,0,0
    if (dstHkey == 'c' ):
        ncol=len( lat_Dst)
    if (dstHkey == 'yx' ):
        ny,nx =len( lat_Dst),len( lon_Dst)

    #if (dstHkey == 'yx') or (srcShape[-2:] == 'yx' ):
    if (srcShape[-2:] == 'yx' ):
        print( "Bomb out ... yx Src shapes not implemented yet " )
        return
    
    print( f" dstHkey={dstHkey} , ncol={ncol} , ny={ny}, nx={nx} " )
    print( f" srcShape={srcShape}, nzot={nzot}, nz={nz} ") 
    
   
    #############################################
    # mesh-to-mesh
    #############################################
    if ((srcShape == 'c' ) and (dstHkey == 'c')) :    
        xfld_Dst = np.zeros( (ncol) , dtype=np.float64 )
        Slice_Src = xfld_Src
        Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                            regrd = regrd , 
                            srcField= srcf , 
                            dstField= dstf , 
                            srcGridkey= srcHkey ,
                            dstGridkey= dstHkey )
        xfld_Dst = Slice_Dst
    if ((srcShape == 'oc' ) and (dstHkey == 'c')) :    
        xfld_Dst = np.zeros( (nzot,ncol) , dtype=np.float64 )
        for tin in np.arange( nzot ):
            Slice_Src = xfld_Src[tin,:]
            Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                regrd = regrd , 
                                srcField= srcf , 
                                dstField= dstf , 
                                srcGridkey= srcHkey ,
                                dstGridkey= dstHkey )
            xfld_Dst[tin,:] = Slice_Dst
    if ((srcShape == 'tzc' ) and (dstHkey == 'c')) :    
        xfld_Dst = np.zeros( (nzot,nz,ncol) , dtype=np.float64 )
        for tin in np.arange( nzot ):
            for zi in np.arange(nz):
                Slice_Src = xfld_Src[tin,zi,:]
                Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                regrd = regrd , 
                                srcField= srcf , 
                                dstField= dstf , 
                                srcGridkey= srcHkey ,
                                dstGridkey= dstHkey )
                xfld_Dst[tin,zi,:] = Slice_Dst

    #############################################
    # mesh-to-grid
    #############################################
    if ((srcShape == 'c' ) and (dstHkey == 'yx')) :    
        xfld_Dst = np.zeros( (ny,nx) , dtype=np.float64 )
        Slice_Src = xfld_Src
        Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                            regrd = regrd , 
                            srcField= srcf , 
                            dstField= dstf , 
                            srcGridkey= srcHkey ,
                            dstGridkey= dstHkey )
        xfld_Dst = Slice_Dst
    if ((srcShape == 'oc' ) and (dstHkey == 'yx')) :    
        xfld_Dst = np.zeros( (nzot,ny,nx) , dtype=np.float64 )
        for tin in np.arange( nzot ):
            Slice_Src = xfld_Src[tin,:]
            Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                regrd = regrd , 
                                srcField= srcf , 
                                dstField= dstf , 
                                srcGridkey= srcHkey ,
                                dstGridkey= dstHkey )
            xfld_Dst[tin,:,:] = Slice_Dst
    if ((srcShape == 'tzc' ) and (dstHkey == 'yx')) :    
        xfld_Dst = np.zeros( (nzot,nz,ny,nx) , dtype=np.float64 )
        for tin in np.arange( nzot ):
            for zi in np.arange(nz):
                Slice_Src = xfld_Src[tin,zi,:]
                Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                regrd = regrd , 
                                srcField= srcf , 
                                dstField= dstf , 
                                srcGridkey= srcHkey ,
                                dstGridkey= dstHkey )
                xfld_Dst[tin,zi,:,:] = Slice_Dst

    return xfld_Dst


############################################################################################################

def Vert(DstVgrid=None, DstTZHkey=None, SrcVgrid=None, xfld_Src=None, ps_Src=None, pmid_output=False ):
    ##############################################################
    #  We will refer to '2D' and '3D' fields regardless of the 
    #  dimensions of the array representing them. '2D' will simply 
    #  mean fields with no vertical or time dimension.
    ############################################################
    DstVgridFile = GrU.gridInfo(Vgrid=DstVgrid,VgridOnly=True) #,Vgrid=DstVgrid)
    SrcVgridFile = GrU.gridInfo(Vgrid=SrcVgrid,VgridOnly=True) #,Vgrid=DstVgrid)

    vDst=xr.open_dataset( DstVgridFile )
    amid_Dst = vDst['hyam'].values
    bmid_Dst = vDst['hybm'].values
    aint_Dst = vDst['hyai'].values
    bint_Dst = vDst['hybi'].values

    vSrc=xr.open_dataset( SrcVgridFile )
    amid_Src = vSrc['hyam'].values
    bmid_Src = vSrc['hybm'].values
    aint_Src = vSrc['hyai'].values
    bint_Src = vSrc['hybi'].values

    p_00 = 100_000. # Here we just use the sensible value of p_00


    pmid_Src,pint_Src,delp_Src \
        = MkP.Pressure (am=amid_Src ,
                        bm=bmid_Src ,
                        ai=aint_Src ,
                        bi=bint_Src ,
                        ps=ps_Src ,
                        p_00=p_00 , 
                        Gridkey = DstTZHkey )

    pmid_Dst,pint_Dst,delp_Dst \
        = MkP.Pressure (am=amid_Dst ,
                        bm=bmid_Dst ,
                        ai=aint_Dst ,
                        bi=bint_Dst ,
                        ps=ps_Src ,
                        p_00=p_00 , 
                        Gridkey = DstTZHkey )

    lnpint_Src = -7_000. * np.log( pint_Src / p_00 )
    lnpmid_Src = -7_000. * np.log( pmid_Src / p_00 )
    lnpint_Dst = -7_000. * np.log( pint_Dst / p_00 )
    lnpmid_Dst = -7_000. * np.log( pmid_Dst / p_00 )

    xfld_Dst = vrg.VertRG( a_x  = xfld_Src ,
                            zSrc = lnpmid_Src ,
                            zDst = lnpmid_Dst ,
                            Gridkey =DstTZHkey ,
                            kind = 'linear' ) #linea
    

    if (pmid_output==True):
        return xfld_Dst,pmid_Dst,pmid_Src
    else:
        return xfld_Dst

