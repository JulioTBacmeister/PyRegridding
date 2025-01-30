#!/usr/bin/env python
import sys
import argparse as arg


workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} ")


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

from PyRegridding.Regridder import scripGen as SG
from PyRegridding.Regridder import esmfRegrid as erg
from PyRegridding.Utils import MyConstants as Con
from PyRegridding.Utils import GridUtils as GrU
from PyRegridding.Utils import MakePressures as MkP
# "ChatGPI version" --- 
from PyRegridding.Regridder import VertRegridFlexLL as vrg
print( "Using Flexible parallel/serial VertRegrid ")

from PyRegridding.Utils import humiditycalcs as hum


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



def Horz(Dst=None,Src=None, xfld_Src=None, RegridMethod=None, RegridObj_In=None, UseFiles=True ):
    
    ##############################################################
    #  We will refer to '2D' and '3D' fields regardless of the 
    #  dimensions of the array representing them. '2D' will simply 
    #  mean fields with no vertical or time dimension.
    ############################################################

    #######################################
    # Grab grid monikers from RegridObj_in 
    # if needed and possible
    #######################################
    if (Src is None) and (Dst is None) and (RegridObj_In is not None):
        if (len(RegridObj_In) >= 5):
            Src=RegridObj_In[3]
            Dst=RegridObj_In[4]
            print( f"Got monikers {Src} and {Dst} from RegridObj_In" )
        else:
            raise ValueError("Need Src and Dst")
    else:
        print( f"Got monikers {Src} and {Dst} from arguments" )

            
    
    DstInfo = GrU.gridInfo(Dst) #,Vgrid=DstVgrid)
    dstHkey = DstInfo['Hkey']
    dst_type =DstInfo['type']
    dst_scrip =DstInfo['scrip']

    SrcInfo = GrU.gridInfo(Src)
    srcHkey = SrcInfo['Hkey']
    src_type =SrcInfo['type']
    src_scrip =SrcInfo['scrip']

        
    
    print( f"\n \n Mapping ... {Src} -x- {Dst} " , flush=True )

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

        regrd, srcf, dstf = erg.GenWrtRdWeights(Dst=Dst , 
                                                Src=Src , 
                                                UseFiles=UseFiles , 
                                                RegridMethod = RegridMethod_  )
        
        if ( xfld_Src is None ):
            print( f" Not interpolating. Returning: regrd, srcf, dstf, {Src}, {Dst} " , flush=True )
            return regrd, srcf, dstf, Src, Dst
            
    else:
        print( f" Getting (regrd, srcf, dstf) from argument " , flush=True ) 
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
        print("your array can't be classified " , flush=True )
        srcShape='nowayjose'
    
    
    ###########################################################

    lat_Dst,lon_Dst = GrU.latlon( scrip=dst_scrip, Hkey=dstHkey )
    print( f"Used GrU.latlon {dst_scrip} for lat lon " , flush=True )

    ncol,ny,nx=0,0,0
    if (dstHkey == 'c' ):
        ncol=len( lat_Dst)
    if (dstHkey == 'yx' ):
        ny,nx =len( lat_Dst),len( lon_Dst)

    #if (dstHkey == 'yx') or (srcShape[-2:] == 'yx' ):
    #if (srcShape[-2:] == 'yx' ) and (dstHkey == 'yx' ):
    #    print( "Bomb out ... yx grid-to-grid not implemented yet " )
    #    return
    
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

    #############################################
    # grid-to-mesh
    #############################################
    if ((srcShape == 'yx' ) and (dstHkey == 'c')) :    
        xfld_Dst = np.zeros( (ncol) , dtype=np.float64 )
        Slice_Src = xfld_Src
        Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                            regrd = regrd , 
                            srcField= srcf , 
                            dstField= dstf , 
                            srcGridkey= srcHkey ,
                            dstGridkey= dstHkey )
        xfld_Dst = Slice_Dst
    if ((srcShape == 'oyx' ) and (dstHkey == 'c')) :    
        xfld_Dst = np.zeros( (nzot,ncol) , dtype=np.float64 )
        for tin in np.arange( nzot ):
            Slice_Src = xfld_Src[tin,:,:]
            Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                regrd = regrd , 
                                srcField= srcf , 
                                dstField= dstf , 
                                srcGridkey= srcHkey ,
                                dstGridkey= dstHkey )
            xfld_Dst[tin,:] = Slice_Dst
    if ((srcShape == 'tzyx' ) and (dstHkey == 'c')) :    
        xfld_Dst = np.zeros( (nzot,nz,ncol) , dtype=np.float64 )
        for tin in np.arange( nzot ):
            for zi in np.arange(nz):
                Slice_Src = xfld_Src[tin,zi,:,:]
                Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                regrd = regrd , 
                                srcField= srcf , 
                                dstField= dstf , 
                                srcGridkey= srcHkey ,
                                dstGridkey= dstHkey )
                xfld_Dst[tin,zi,:] = Slice_Dst

    #############################################
    # grid-to-grid
    #############################################
    if ((srcShape == 'yx' ) and (dstHkey == 'yx')) :    
        xfld_Dst = np.zeros( (ny,nx) , dtype=np.float64 )
        Slice_Src = xfld_Src
        Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                            regrd = regrd , 
                            srcField= srcf , 
                            dstField= dstf , 
                            srcGridkey= srcHkey ,
                            dstGridkey= dstHkey )
        xfld_Dst = Slice_Dst
    if ((srcShape == 'oyx' ) and (dstHkey == 'yx')) :    
        xfld_Dst = np.zeros( (nzot,ny,nx) , dtype=np.float64 )
        for tin in np.arange( nzot ):
            Slice_Src = xfld_Src[tin,:,:]
            Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                regrd = regrd , 
                                srcField= srcf , 
                                dstField= dstf , 
                                srcGridkey= srcHkey ,
                                dstGridkey= dstHkey )
            xfld_Dst[tin,:,:] = Slice_Dst
    if ((srcShape == 'tzyx' ) and (dstHkey == 'yx')) :    
        xfld_Dst = np.zeros( (nzot,nz,ny,nx) , dtype=np.float64 )
        for tin in np.arange( nzot ):
            for zi in np.arange(nz):
                Slice_Src = xfld_Src[tin,zi,:,:]
                Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                                regrd = regrd , 
                                srcField= srcf , 
                                dstField= dstf , 
                                srcGridkey= srcHkey ,
                                dstGridkey= dstHkey )
                xfld_Dst[tin,zi,:,:] = Slice_Dst

    return xfld_Dst


############################################################################################################
def HorzSlice(Dst=None, Src=None, xfld_Src=None, RegridMethod=None, RegridObj_In=None ):
    ##################################################################################
    # Bare-bones horizontal regridding of ONE horizontal (physically speaking)
    # slice of data
    ##################################################################################

    # If ESMF regridding objects are not provided, then generate them.
    # Note, in this case Src and Dst 'nicknames' for horizontal grids
    # must also be provided.
    #--------------------
    if ( RegridObj_In == None ):
        
        # Get Source grid Info
        SrcInfo = GrU.gridInfo(Src)
        srcHkey = SrcInfo['Hkey']
        src_type =SrcInfo['type']
        src_scrip =SrcInfo['scrip']
        
        # Get Destination grid Info
        DstInfo = GrU.gridInfo(Dst) 
        dstHkey = DstInfo['Hkey']
        dst_type =DstInfo['type']
        dst_scrip =DstInfo['scrip']
        
        # Decide on method
        if (RegridMethod==None):
            RegridMethod_ = 'CONSERVE_2ND' 
        else:
            RegridMethod_ = RegridMethod
            
        # Generate ESMF objects
        regrd, srcf, dstf = erg.GenWrtRdWeights(Dst=Dst , 
                                                Src=Src , 
                                                UseFiles=True , 
                                                RegridMethod = RegridMethod_  )
        # If condition==True return ESMF objects in a tuple (exit)       
        if ( xfld_Src is None ):
            print( f" Not interpolating. Returning: regrd, srcf, dstf ", flush=True  )
            return regrd, srcf, dstf
            
    # ESMF regridding objects have been provided in a tuple.
    #--------------------
    else:
        print( f" Getting (regrd, srcf, dstf) from argument ", flush=True ) 
        regrd, srcf, dstf = RegridObj_In[0], RegridObj_In[1], RegridObj_In[2]   

    
    #############################################
    # mesh-to-mesh
    #############################################
    if ((srcHkey == 'c' ) and (dstHkey == 'c')) :    
        srcf.data = xfld_Src
        r  = regrd(  srcf ,  dstf )
        xfld_Dst = dstf.data 
    #############################################
    # mesh-to-grid
    #############################################
    if ((srcHkey == 'c' ) and (dstHkey == 'xy')) :    
        srcf.data = xfld_Src
        r  = regrd(  srcf ,  dstf )
        xfld_Dst = dstf.data.transpose() 
    #############################################
    # grid-to-mesh
    #############################################
    if ((srcHkey == 'xy' ) and (dstHkey == 'c')) :    
        srcf.data = xfld_Src.transpose()
        r  = regrd(  srcf ,  dstf )
        xfld_Dst = dstf.data 
    #############################################
    # grid-to-grid
    #############################################
    if ((srcHkey == 'xy' ) and (dstHkey == 'c')) :    
        srcf.data = xfld_Src.transpose()
        r  = regrd(  srcf ,  dstf )
        xfld_Dst = dstf.data.transpose()
        

    return xfld_Dst


############################################################################################################

def Vert(DstVgrid=None, DstTZHkey=None, SrcVgrid=None, xfld_Src=None, ps_Src=None, pmid_output=False ):

    if ( (DstVgrid == SrcVgrid) and (xfld_Src is not None) ):
        print( f" No need to interpolate in the vertical " )
        xfld_Dst = xfld_Src
        return xfld_Dst

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

