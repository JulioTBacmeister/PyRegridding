# Import packages 
import sys
workdir_ = '/glade/work/juliob'
if ( workdir_ not in sys.path ):
    sys.path.append(workdir_)
    print( f" a path to {workdir_} added in {__name__} ")

import xarray as xr
import numpy as np

"""
try:
    import ESMF as E
except ImportError:
    import esmpy as E
"""

import esmpy as E
import importlib

import copy
import glob


def Regrid( srcScrip , dstScrip , srcType , dstType ,  RegridMethod="CONSERVE" , **kwargs ):

    write_wgts = False
    read_wgts  = False
    if ("write_weights" in kwargs):
        write_wgts=kwargs["write_weights"]
        wgts_file =kwargs["weights_file"]
    if ("read_weights" in kwargs):
        read_wgts=kwargs["read_weights"]
        wgts_file =kwargs["weights_file"]

    if(RegridMethod.upper()=='CONSERVE'):
        regrid_method=E.RegridMethod.CONSERVE 
    if(RegridMethod.upper()=='CONSERVE_2ND'):
        regrid_method=E.RegridMethod.CONSERVE_2ND 
    if(RegridMethod.upper()=='BILINEAR'):
        regrid_method=E.RegridMethod.BILINEAR 
    
    if(srcType.lower()=='mesh'):
        srcDesc=E.Mesh( filename=srcScrip,
            filetype=E.FileFormat.SCRIP )
        
    if(srcType.lower()=='grid'):
        srcDesc=E.Grid( filename=srcScrip ,
            filetype=E.FileFormat.SCRIP ,
            add_corner_stagger=True   )
        
    if(dstType.lower()=='mesh'):
        dstDesc=E.Mesh( filename=dstScrip,
            filetype=E.FileFormat.SCRIP )
        
    if(dstType.lower()=='grid'):
        dstDesc=E.Grid( filename=dstScrip ,
            filetype=E.FileFormat.SCRIP ,
            add_corner_stagger=True   )
        
        
    if(srcType.lower()=='mesh'):
        srcField = E.Field(srcDesc, meshloc=E.MeshLoc.ELEMENT )
    if(srcType.lower()=='grid'):
        srcField = E.Field(srcDesc )

    if(dstType.lower()=='mesh'):
        dstField = E.Field(dstDesc, meshloc=E.MeshLoc.ELEMENT )
    if(dstType.lower()=='grid'):
        dstField = E.Field(dstDesc )

    
    dstField.data[:]=1e20
    srcField.data[:]=1e20

    if (write_wgts==True):
        print(f"Generating regridding weights. Method {RegridMethod} : ESMF method= {regrid_method}")
        print(f"And writing to {wgts_file}")
        Regrd = E.Regrid( srcField , dstField , 
                      filename = wgts_file,
                      regrid_method=regrid_method,
                      unmapped_action=E.UnmappedAction.IGNORE)
    elif (read_wgts==True):
        print(f"Reading weights file {wgts_file}")
        Regrd = E.RegridFromFile( srcField , dstField , 
                                 filename=wgts_file )
    else:
        print(f"Generating regridding weights. Method {RegridMethod} : ESMF method= {regrid_method}")
        Regrd = E.Regrid( srcField , dstField , 
                      regrid_method=regrid_method,
                      unmapped_action=E.UnmappedAction.IGNORE)


    return Regrd, srcField , dstField 

#########################
#def HorzRG( aSrc, regrd , srcField , dstField , srcShape, dstShape, srcGridkey, dstGridkey ):
def HorzRG( aSrc, regrd , srcField , dstField , srcGridkey, dstGridkey ):

    import copy
    
    #---------------------------------------------------------------------
    # This function takes input ndarray aSrc and remaps in the HORIZONTAL
    # to the destination grid.  Input aSrc must contain at least one 
    # horizontal slice (Hslice), but can also be shaped 
    # Z x Hslice or T x Z x Hslice where Z and T refer to vertical
    # time dimensions.
    # regrd, srcField and dstField are precomupted regridding objects
    # {src,dst}Grid are strings like 'tzyx',where:
    #          t -> time
    #          z -> vertical
    #-----------------------------------------------------------------------
    #------------------------------------------------------------------------
    # We shouldn't really need srcShape and dstShape arguments as long as
    # as have the 'Gridkeys'. The shapes should avaiable from
    #
    #      srcShape = np.shape( aSrc ) 
    #      dstShape = np.shape( dstField.data ) {transpose if 'yx'} 
    #------------------------------------------------------------------------
    
    #srcShape = np.shape( srcField.data )
    srcShape = np.shape( aSrc )
    dstShape = np.shape( dstField.data )
    
    
    #------------------------------------------------------------------------
    # Following block deals with logically-rectangular 'yx' input Horz grids
    #------------------------------------------------------------------------
    if (srcGridkey == 'yx' ):
        srcField.data[:,:] = aSrc.transpose()
        if (dstGridkey=='c'):
            r  = regrd(  srcField,  dstField )
            aDst = copy.deepcopy( dstField.data[:] )
        if (dstGridkey=='yx'):
            r  = regrd(  srcField,  dstField )
            aDst = copy.deepcopy( dstField.data[:,:].transpose() )
            
    if (srcGridkey == 'zyx' ):
        nlev = srcShape[0]
        if (dstGridkey=='c'):
            ncol = dstShape[0]
            aDst = np.zeros([nlev,ncol])
            for L in np.arange(nlev):
                srcField.data[:,:] = aSrc[L,:,:].transpose()
                r  = regrd(  srcField,  dstField )
                aDst[L,:] = copy.deepcopy( dstField.data[:] )
        if (dstGridkey=='yx'):
            ny,nx = dstShape[1],dstShape[0]
            aDst = np.zeros([nlev,ny,nx])
            for L in np.arange(nlev):
                srcField.data[:,:] = aSrc[L,:,:].transpose()
                r  = regrd(  srcField,  dstField )
                aDst[L,:,:] = copy.deepcopy( dstField.data[:,:].transpose() )
                
    if (srcGridkey == 'tyx' ):
        ntim = srcShape[0]
        if (dstGridkey=='c'):
            ncol = dstShape[0]
            aDst = np.zeros([ntim,ncol])
            for i in np.arange(ntim):
                srcField.data[:,:] = aSrc[i,:,:].transpose()
                r  = regrd(  srcField,  dstField )
                aDst[i,:] = copy.deepcopy( dstField.data[:] )
        if (dstGridkey=='yx'):
            ny,nx = dstShape[1],dstShape[0]
            aDst = np.zeros([ntim,ny,nx])
            for i in np.arange(ntim):
                srcField.data[:,:] = aSrc[i,:,:].transpose()
                r  = regrd(  srcField,  dstField )
                aDst[i,:,:] = copy.deepcopy( dstField.data[:,:].transpose() )
                
    if (srcGridkey == 'tzyx' ):
        ntim = srcShape[0]
        nlev = srcShape[1]
        if (dstGridkey=='c'):
            ncol = dstShape[0]
            aDst = np.zeros([ntim,nlev,ncol])
            for i in np.arange(ntim):
                for L in np.arange(nlev):
                    srcField.data[:,:] = aSrc[i,L,:,:].transpose()
                    r  = regrd(  srcField,  dstField )
                    aDst[i,L,:] = copy.deepcopy( dstField.data[:] )
        if (dstGridkey=='yx'):
            ny,nx = dstShape[1],dstShape[0]
            aDst = np.zeros([ntim,nlev,ny,nx])
            for i in np.arange(ntim):
                for L in np.arange(nlev):
                    srcField.data[:,:] = aSrc[i,L,:,:].transpose()
                    r  = regrd(  srcField,  dstField )
                    aDst[i,L,:,:] = copy.deepcopy( dstField.data[:,:].transpose() )
                    
    #-----------------------------------------------------------------
    # Following block deals with unstructured 'c' input Horz grids
    #------------------------------------------------------------------
    if (srcGridkey == 'c' ):
        srcField.data[:] = aSrc[:]
        if (dstGridkey=='c'):
            r  = regrd(  srcField,  dstField )
            aDst = copy.deepcopy( dstField.data[:] )
        if (dstGridkey=='yx'):
            r  = regrd(  srcField,  dstField )
            aDst = copy.deepcopy( dstField.data[:,:].transpose() )
            
    if (srcGridkey == 'zc' ):
        nlev = srcShape[0]
        if (dstGridkey=='c'):
            ncol = dstShape[0]
            aDst = np.zeros([nlev,ncol])
            for L in np.arange(nlev):
                srcField.data[:] = aSrc[L,:]
                r  = regrd(  srcField,  dstField )
                aDst[L,:] = copy.deepcopy( dstField.data[:] )
        if (dstGridkey=='yx'):
            ny,nx = dstShape[1],dstShape[0]
            aDst = np.zeros([nlev,ny,nx])
            for L in np.arange(nlev):
                srcField.data[:] = aSrc[L,:]
                r  = regrd(  srcField,  dstField )
                aDst[L,:,:] = copy.deepcopy( dstField.data[:,:].transpose() )
                
    if (srcGridkey == 'tc' ):
        ntim = srcShape[0]
        if (dstGridkey=='c'):
            ncol = dstShape[0]
            aDst = np.zeros([ntim,ncol])
            for i in np.arange(ntim):
                srcField.data[:] = aSrc[i,:]
                r  = regrd(  srcField,  dstField )
                aDst[i,:] = copy.deepcopy( dstField.data[:] )
        if (dstGridkey=='yx'):
            ny,nx = dstShape[1],dstShape[0]
            aDst = np.zeros([ntim,ny,nx])
            for i in np.arange(ntim):
                srcField.data[:,:] = aSrc[i,:]
                r  = regrd(  srcField,  dstField )
                aDst[i,:,:] = copy.deepcopy( dstField.data[:,:].transpose() )
                
    if (srcGridkey == 'tzc' ):
        ntim = srcShape[0]
        nlev = srcShape[1]
        if (dstGridkey=='c'):
            ncol = dstShape[0]
            aDst = np.zeros([ntim,nlev,ncol])
            for i in np.arange(ntim):
                for L in np.arange(nlev):
                    srcField.data[:] = aSrc[i,L,:]
                    r  = regrd(  srcField,  dstField )
                    aDst[i,L,:] = copy.deepcopy( dstField.data[:] )
        if (dstGridkey=='yx'):
            ny,nx = dstShape[1],dstShape[0]
            aDst = np.zeros([ntim,nlev,ny,nx])
            for i in np.arange(ntim):
                for L in np.arange(nlev):
                    srcField.data[:] = aSrc[i,L,:]
                    r  = regrd(  srcField,  dstField )
                    aDst[i,L,:,:] = copy.deepcopy( dstField.data[:,:].transpose() )
                    
    return aDst

######################################################################
#def WriteWeights( srcScrip , dstScrip , srcType , dstType ,  RegridMethod="CONSERVE_2ND" , **kwargs ):
def GenWrtRdWeights( Src=None , Dst=None ,  RegridMethod="CONSERVE_2ND" , UseFiles=None, **kwargs ):
    """
    This function returns the regridding objects (Regrid, srcField, dstField) used by ESMF.
    Src, Dst inputs (required) are code names for grids as recognized in 
        PyRegridding/Utils/GridUtils.py [.gridInfo]
    UseFiles (required) decides whether weight files will be processed. 3 behaviors possible.
        UseFiles=False
            Simply generate the emsf regridding objects and return
        UseFiles=True
            wgts_file string is generated from Dst,Src, RegridMethod
            If a file name=wgts_file exists
                Read weights and return esmf regridding objects
            If a file name=wgts_file DOES NOT exist
                Generate weights, write file and return esmf regridding objects
    RegridMethod (optional)
    """

    import time
    import os
    from PyRegridding.Utils import GridUtils as GrU

    ##############################
    #  Begin ...
    ##############################

    
    tic_overall = time.perf_counter()

    if (UseFiles is None ):
        ErrString = f'You ABSOLUTELY have to set UseFiles = True or False'
        print( ErrString )
        return ErrString,-99,-99

        
    DstInfo  = GrU.gridInfo(Dst) #,Vgrid=DstVgrid)
    dstHkey  = DstInfo['Hkey']
    dstType  =DstInfo['type']
    dstScrip =DstInfo['scrip']

    SrcInfo = GrU.gridInfo(Src)
    srcHkey = SrcInfo['Hkey']
    srcType =SrcInfo['type']
    srcScrip =SrcInfo['scrip']
    
    if(RegridMethod.upper()=='CONSERVE'):
        regrid_method=E.RegridMethod.CONSERVE 
    if(RegridMethod.upper()=='CONSERVE_2ND'):
        regrid_method=E.RegridMethod.CONSERVE_2ND 
    if(RegridMethod.upper()=='BILINEAR'):
        regrid_method=E.RegridMethod.BILINEAR 

    toc_here = time.perf_counter()
    pTime = f"Ready to go ... {Src} -x- {Dst} {toc_here - tic_overall:0.4f} seconds"
    print(pTime, flush=True )
    
    if(srcType.lower()=='mesh'):
        srcDesc=E.Mesh( filename=srcScrip,
            filetype=E.FileFormat.SCRIP )
        
    if(srcType.lower()=='grid'):
        srcDesc=E.Grid( filename=srcScrip ,
            filetype=E.FileFormat.SCRIP ,
            add_corner_stagger=True   )
        
    if(dstType.lower()=='mesh'):
        dstDesc=E.Mesh( filename=dstScrip,
            filetype=E.FileFormat.SCRIP )
        
    if(dstType.lower()=='grid'):
        dstDesc=E.Grid( filename=dstScrip ,
            filetype=E.FileFormat.SCRIP ,
            add_corner_stagger=True   )
        
    toc_here = time.perf_counter()
    pTime = f"Finished Grid Mesh creation  {toc_here - tic_overall:0.4f} seconds"
    print(pTime, flush=True )
        
        
    if(srcType.lower()=='mesh'):
        srcField = E.Field(srcDesc, meshloc=E.MeshLoc.ELEMENT )
    if(srcType.lower()=='grid'):
        srcField = E.Field(srcDesc )

    if(dstType.lower()=='mesh'):
        dstField = E.Field(dstDesc, meshloc=E.MeshLoc.ELEMENT )
    if(dstType.lower()=='grid'):
        dstField = E.Field(dstDesc )

    print( dstScrip ,'\n',np.shape(dstField.data) )   
    dstField.data[:]=1e20
    print( srcScrip ,'\n',np.shape(srcField.data) )   
    srcField.data[:]=1e20

    toc_here = time.perf_counter()
    pTime = f"Finished dstField,srcField  {toc_here - tic_overall:0.4f} seconds"
    print(pTime, flush=True )

    if (UseFiles==True):
        wgts_file=f'/glade/work/juliob/GridFiles/Weights/{Src}_x_{Dst}_{RegridMethod.upper()}.nc'
        print( wgts_file)
        if (not os.path.exists(wgts_file)):
            print(f"Generating regridding weights. Method {RegridMethod} : ESMF method= {regrid_method}")
            Regrd = E.Regrid( srcField , dstField , 
                          filename = wgts_file,
                          regrid_method=regrid_method,
                          unmapped_action=E.UnmappedAction.IGNORE)
            toc_here = time.perf_counter()
            pTime = f"Finished generatind Regrd from scratch {toc_here - tic_overall:0.4f} seconds"
            print(pTime , flush=True )
        else:
            print(f"Reading weights from {wgts_file} ")
            Regrd = E.RegridFromFile( srcField , dstField , 
                          filename = wgts_file )
            toc_here = time.perf_counter()
            pTime = f"Finished generatind Regrd from wgts file {toc_here - tic_overall:0.4f} seconds"
            print(pTime , flush=True )
    else:
        print(f"Not dealing with weight files at all - just calculating Regrid object")
        Regrd = E.Regrid( srcField , dstField , 
                      regrid_method=regrid_method,
                      unmapped_action=E.UnmappedAction.IGNORE)
        
        toc_here = time.perf_counter()
        pTime = f"Finished generatind Regrd from scratch {toc_here - tic_overall:0.4f} seconds"
        print(pTime , flush=True )
      
    return Regrd, srcField , dstField 

#########################

