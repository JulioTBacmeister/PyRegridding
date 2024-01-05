#!/usr/bin/env python
# Import packages 

import sys
sys.path.append('/glade/work/juliob/PyRegridding/Utils/')

import numpy as np

import GridUtils as GU
import esmfRegrid as erg
import latlon_w_scrip as LL

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
    
    src_scrip, src_Hkey, src_type = GU.scrip_etc( grid=agrid )
    dst_scrip, dst_Hkey, dst_type = GU.scrip_etc( grid=bgrid )
    
    lat_dst,lon_dst = LL.latlon( scrip=dst_scrip , gridHkey = dst_Hkey )

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
