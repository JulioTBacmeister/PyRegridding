#!/usr/bin/env python
# Import packages 

import sys
import os
#sys.path.append('/glade/work/juliob/PyRegridding/Utils/')

#-----------------------------------------
# Find path to this module and calc/append 
# paths relative to this path 
#------------------------------------------
module_a_dir = os.path.dirname(os.path.abspath(__file__))
utils_path = os.path.join(module_a_dir, '..', 'Utils')
sys.path.append(utils_path)
print( f" a path added in {__name__} {utils_path} ")


import numpy as np
from scipy import interpolate as intr

import GridUtils as GU
import esmfRegrid as erg
import time

import importlib
importlib.reload( GU )


def Hregrid(avar=None,
            agrid=None,
            akey=None,
            bgrid=None,
            RegridMethod='BILINEAR',
            CreateRegrid=False,
            regrd=None,
            dstf=None,
            srcf=None):
    ####################
    # avar (float): array to regridded
    # agrid(char): desc code for grid avar lives on, e.g., 'ne30pg3' 
    # akey(char): shape of input variable to bre horz regridded
    #        Possible values 'c','yx','zc','zyx','tzc',tzyx'
    # bgrid(char): code for taarget grid. e.g., 'fv0.9x1.25'
    #  
    ####################
    
    #src_scrip, src_Hkey, src_type = GU.gridInfo( grid=agrid )
    #dst_scrip, dst_Hkey, dst_type = GU.gridInfo( grid=bgrid )
    AgridInfo = GU.gridInfo( grid=agrid )
    src_scrip, src_Hkey, src_type = AgridInfo['scrip'],AgridInfo['Hkey'],AgridInfo['type']
    
    BgridInfo = GU.gridInfo( grid=bgrid )
    dst_scrip, dst_Hkey, dst_type = BgridInfo['scrip'],BgridInfo['Hkey'],BgridInfo['type']

    print(f" Agrid {agrid} {src_scrip} , {src_Hkey} , {src_type}")
    print(f" Bgrid {bgrid} {dst_scrip} , {dst_Hkey} , {dst_type}")

    
    lat_dst,lon_dst = GU.latlon( scrip=dst_scrip , Hkey = dst_Hkey )

    if ( dst_Hkey == 'c' ):
        nc_dst = len( lat_dst )
    if ( dst_Hkey == 'yx' ):
        ny_dst, nx_dst = len( lat_dst ) , len( lon_dst )
    

    if ( regrd == None ):
        # ----------------------------------------------
        # Make object for ESMF regridding from SRC
        # grid to CAM target. Scrip files need to be provided even 
        # when a weight file is used.
        # We make both conservative and bilinear mapping
        # files, because we will use conservative for 2D surface vars, but
        # bilinear for 3D met fields.
        # ----------------------------------------------
        regrd, srcf, dstf = erg.Regrid( srcScrip = src_scrip , 
                                        srcType  = src_type  ,
                                        dstScrip = dst_scrip ,
                                        dstType  = dst_type  ,
                                        RegridMethod = RegridMethod )
    else:
        print( 'Need to have regrd,srcf, and dstf as args ' )

    if ( CreateRegrid == True ):
        return regrd,srcf,dstf

    ####################
    # The following 2D vs 3D determination
    # should be done in a better way.
    # What is here will work for SE CAM output, i.e., 
    # (time,col) or (time,lev,col)
    #######################
    # Remapping of purely 2D (Horz) vars
    if ( akey == 'c' or akey == 'yx' ):    
        #xfld_Dst = np.zeros( (nt,ny,nx) , dtype=np.float64 )
        nlev=1
        Slice_Src = avar
        Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                     regrd = regrd , 
                     srcField= srcf , 
                     dstField= dstf , 
                     srcGridkey= src_Hkey ,
                     dstGridkey= dst_Hkey )
        bvar = Slice_Dst

    #############
    # Remapping for 3D vars {height,time} x Horz
    if ( akey == 'tc' or akey == 'tyx' or akey == 'zc' or akey == 'zyx' ):    
        nlev = np.shape( avar )[0] # no of slices 
        if (dst_Hkey == 'c'):
            bvar = np.zeros( (nlev , nc_dst) , dtype=np.float64 )
        if (dst_Hkey == 'yx'):
            bvar = np.zeros( (nlev , ny_dst, nx_dst ) , dtype=np.float64 )
            
        for L in np.arange( nlev ):
            if (src_Hkey == 'c'):
                Slice_Src = avar[L,:]
            if (src_Hkey == 'yx'):
                Slice_Src = avar[L,:,:]
            
            Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                         regrd = regrd , 
                         srcField= srcf , 
                         dstField= dstf , 
                         srcGridkey= src_Hkey ,
                         dstGridkey= dst_Hkey )
            
            if (dst_Hkey == 'c'):
                bvar[L,:] = Slice_Dst
            if (dst_Hkey == 'yx'):
                bvar[L,:,:] = Slice_Dst

    #############
    # Remapping for 4D vars time x height x Horz
    if ( akey == 'tzc' or akey == 'tzyx'):    
        ntim,nlev = np.shape( avar )[0],np.shape( avar )[1] # no of slices 
        if (dst_Hkey == 'c'):
            bvar = np.zeros( (ntim, nlev , nc_dst) , dtype=np.float64 )
        if (dst_Hkey == 'yx'):
            bvar = np.zeros( (ntim, nlev , ny_dst, nx_dst ) , dtype=np.float64 )
            
        for n in np.arange( ntim ):
            for L in np.arange( nlev ):
                if (src_Hkey == 'c'):
                    Slice_Src = avar[n,L,:]
                if (src_Hkey == 'yx'):
                    Slice_Src = avar[n,L,:,:]

                Slice_Dst = erg.HorzRG( aSrc = Slice_Src , 
                             regrd = regrd , 
                             srcField= srcf , 
                             dstField= dstf , 
                             srcGridkey= src_Hkey ,
                             dstGridkey= dst_Hkey )

                if (dst_Hkey == 'c'):
                    bvar[n,L,:] = Slice_Dst
                if (dst_Hkey == 'yx'):
                    bvar[n,L,:,:] = Slice_Dst


        
    return bvar,lat_dst,lon_dst

def interpolate_column(zSrcT_col, a_xT_col, zDstT_col, fill_value, kind):
    """Interpolate a single column of data."""
    fint = intr.interp1d(x=zSrcT_col, y=a_xT_col, fill_value=fill_value, kind=kind)
    return fint(zDstT_col)


def VertRG( a_x , zSrc, zDst, Gridkey , fill_value='extrapolate', kind='linear' ):

    # Initialize performance counter
    tic = time.perf_counter()


    #--------------------------------------------------------------------------------
    # Assumes shapes of a_x, zSrc are the same and conformable with the shape of zDst.
    # That is, the horizontal shape of a_x,zSrc and zDst are the same.
    #--------------------------------------------------------------------------------
    # It seems like a better idea to do all the reshping on entry and just maintain 
    # one regridding loop. So, at some point implement reshaping to 'tzc' at top.
    #--------------------------------------------------------------------------------
    
    if (Gridkey == 'zc'):
        nzS,ncol = np.shape( zSrc )
        nzD,ncol = np.shape( zDst )
        a_xz = np.zeros( (nzD,ncol) )
        # Serial loop (Faster if nworkers = 1)
        for i in np.arange(ncol):
            fint=intr.interp1d( x = zSrc[:,i], y=a_x[:,i] , 
                                fill_value=fill_value, kind=kind  )
            a_xz[:, i] = fint(   zDst[:, i ] )
                                  

    if (Gridkey == 'tzc'):
        # Reshape arrays
        nt,nzS,ncol = np.shape( zSrc )
        nt,nzD,ncol = np.shape( zDst )
        a_xT        = np.reshape( np.transpose( a_x , (1,0,2) ) , (nzS,nt*ncol) )
        zSrcT       = np.reshape( np.transpose( zSrc , (1,0,2) ) , (nzS,nt*ncol) )
        zDstT       = np.reshape( np.transpose( zDst , (1,0,2) ) , (nzD,nt*ncol) )
        nzS,ntcol = np.shape( zSrcT )
        nzD,ntcol = np.shape( zDstT )
        a_xzT = np.zeros( (nzD,ntcol) )
        
        # Serial loop (Faster if nworkers = 1)
        for i in np.arange(ntcol):
            fint=intr.interp1d( x = zSrcT[:,i], y=a_xT[:,i] , 
                                fill_value=fill_value, kind=kind  )
            a_xzT[:, i] = fint(   zDstT[:, i ] )


        a_xz = np.transpose( np.reshape( a_xzT, (nzD,nt,ncol) ), (1,0,2) )

    if (Gridkey == 'tzyx'):
        # Reshape arrays
        nt,nzS,ny,nx = np.shape( zSrc )
        nt,nzD,ny,nx = np.shape( zDst )
        a_xT        = np.reshape( np.transpose( a_x , (1,0,2,3) ) , (nzS,nt*nx*ny) )
        zSrcT       = np.reshape( np.transpose( zSrc , (1,0,2,3) ) , (nzS,nt*nx*ny) )
        zDstT       = np.reshape( np.transpose( zDst , (1,0,2,3) ) , (nzD,nt*nx*ny) )
        nzS,ntcol = np.shape( zSrcT )
        nzD,ntcol = np.shape( zDstT )
        a_xzT = np.zeros( (nzD,ntcol) )
        
        # Serial loop (Faster if nworkers = 1)
        for i in np.arange(ntcol):
            fint=intr.interp1d( x = zSrcT[:,i], y=a_xT[:,i] , 
                                fill_value=fill_value, kind=kind  )
            a_xzT[:, i] = fint(   zDstT[:, i ] )

        a_xz = np.transpose( np.reshape( a_xzT, (nzD,nt,ny,nx) ), (1,0,2,3) )

    if (Gridkey == 'zyx'):
        # Reshape arrays
        nzS,ny,nx = np.shape( zSrc )
        nzD,ny,nx = np.shape( zDst )
        a_xT        = np.reshape( a_x , (nzS, nx*ny) )
        zSrcT       = np.reshape( zSrc , (nzS, nx*ny) )
        zDstT       = np.reshape( zDst  , (nzD, nx*ny) )
        nzS,ntcol = np.shape( zSrcT )
        nzD,ntcol = np.shape( zDstT )
        a_xzT = np.zeros( (nzD,ntcol) )
        
        # Serial loop (Faster if nworkers = 1)
        for i in np.arange(ntcol):
            fint=intr.interp1d( x = zSrcT[:,i], y=a_xT[:,i] , 
                                fill_value=fill_value, kind=kind  )
            a_xzT[:, i] = fint(   zDstT[:, i ] )

        a_xz = np.reshape( a_xzT, (nzD,ny,nx) ) 
        
    toc = time.perf_counter()
    IntrTime = f"Pll'zd Vertical int {toc - tic:0.4f} seconds"
    print(IntrTime)
        
        
    return a_xz

