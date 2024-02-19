# Import packages 
import sys
import os
#sys.path.append('../Utils/')
#sys.path.append('/glade/work/juliob/PyRegridding/Utils/')
#-----------------------------------------
# Find path to this module and calc/append 
# paths relative to this path 
#------------------------------------------
module_a_dir = os.path.dirname(os.path.abspath(__file__))
utils_path = os.path.join(module_a_dir, '..', 'Utils')
sys.path.append(utils_path)
print( f" a path added in {__name__} {utils_path} ")


import xarray as xr
import numpy as np
from scipy import interpolate as intr


try:
    import ESMF as E
except ImportError:
    import esmpy as E

import importlib
import glob
import copy
import time

import scripGen as SG
import esmfRegrid as erg
import MyConstants as Con


importlib.reload( erg )
importlib.reload(SG)

import multiprocessing 
#from multiprocessing import Pool, cpu_count  #()

"""
class scipy.interpolate.interp1d(x, y, kind='linear', axis=-1, copy=True, bounds_error=None, fill_value=nan, assume_sorted=False)[source]
"""

Rdry = Con.Rdry()
grav = Con.grav()

def interpolate_column(zSrcT_col, a_xT_col, zDstT_col, fill_value, kind):
    """Interpolate a single column of data."""
    fint = intr.interp1d(x=zSrcT_col, y=a_xT_col, fill_value=fill_value, kind=kind)
    return fint(zDstT_col)


def VertRG( a_x , zSrc, zDst, Gridkey , fill_value='extrapolate', kind='linear' ,npools=1):

    # Initialize performance counter
    tic = time.perf_counter()

    # set number of pools to create
    #nworkers = cpu_count()
    nworkers = len(os.sched_getaffinity(0))
    print('Flexible VertRegrid using reshaped ARRAYS ')
    print('nworkers available ',nworkers)
    if (nworkers>1):
        print('full loop pll-ized according to ChatGPT')
    else:
        print('Using serial code which is faster for 1 worker')
    # setting npools to nworkers might not be considerate of others...
    npools = nworkers


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
        if (nworkers > 1):
            # parallelized loop
            """Interpolate all columns of data in parallel."""
            with multiprocessing.Pool(nworkers) as pool1:
                results = pool1.starmap(
                    interpolate_column,
                    [(zSrc[:, i], a_x[:, i], zDst[:, i], fill_value, kind)
                     for i in range(zSrc.shape[1])] )
                a_xz = np.column_stack(results) 
        else:
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
        
        if (nworkers > 1):
            # parallelized loop
            # Interpolate all columns of data in parallel
            with multiprocessing.Pool(nworkers) as pool1:
                results = pool1.starmap(
                    interpolate_column,
                    [(zSrcT[:, i], a_xT[:, i], zDstT[:, i], fill_value, kind)
                     for i in range(zSrcT.shape[1])] )
                a_xzT = np.column_stack(results)
        else:
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
        
        if (nworkers > 1):
            # parallelized loop
            # Interpolate all columns of data in parallel
            with multiprocessing.Pool(nworkers) as pool1:
                results = pool1.starmap(
                    interpolate_column,
                    [(zSrcT[:, i], a_xT[:, i], zDstT[:, i], fill_value, kind)
                     for i in range(zSrcT.shape[1])]  )
                a_xzT = np.column_stack(results) 
        else:
            # Serial loop (Faster if nworkers = 1)
            for i in np.arange(ntcol):
                fint=intr.interp1d( x = zSrcT[:,i], y=a_xT[:,i] , 
                                    fill_value=fill_value, kind=kind  )
                a_xzT[:, i] = fint(   zDstT[:, i ] )

        a_xz = np.transpose( np.reshape( a_xzT, (nzD,nt,ny,nx) ), (1,0,2,3) )
        
    toc = time.perf_counter()
    IntrTime = f"Pll'zd Vertical int {toc - tic:0.4f} seconds"
    print(IntrTime)
        
        
    return a_xz


def TsExtrap( ps, pm150, te150 ):

    dTdz_std = -6.5 * 1.e-3

    ts_x = te150 - dTdz_std * (Rdry/grav) * ( ( ps/pm150 ) - 1. ) * te150

    return ts_x

def TeWO (te ,pmid, te150, pm150, ts, ps, L150, Gridkey ):

    teWO = copy.deepcopy(te)
    lnpmid = np.log( pmid )
    lnps = np.log( ps )
    lnpm150 = np.log( pm150 )

    if ( Gridkey == 'tzc' ):
        nt,nz,ncol = np.shape( teWO )
        for i in np.arange( nt ):
            for c in np.arange( ncol ):
                for z in np.arange( start=L150[i,c], stop=nz, step=1 ):
                    dtdlnp = ( te150[i,c] - ts[i,c] ) / ( lnpm150[i,c] - lnps[i,c] )
                    teWO[i,z,c] = te150[i,c] + dtdlnp * ( lnpmid[i,z,c] - lnpm150[i,c] )
                        
    if ( Gridkey == 'tzyx' ):
        nt,nz,ny,nx = np.shape( teWO )
        for i in np.arange( nt ):
            for y in np.arange( ny ):
                for x in np.arange( nx ):
                    for z in np.arange( start=L150[i,y,x], stop=nz, step=1 ):
                        dtdlnp = ( te150[i,y,x] - ts[i,y,x] ) / ( lnpm150[i,y,x] - lnps[i,y,x] )
                        teWO[i,z,y,x] = te150[i,y,x] + dtdlnp * ( lnpmid[i,z,y,x] - lnpm150[i,y,x] )
                        
        

    return teWO

def BottomFill (a_zCAM , a_zERA, pmid_zCAM, ps_ERA, Gridkey ):

    tic = time.perf_counter()

    a_zCAMf = copy.deepcopy(a_zCAM)

    if ( Gridkey == 'tzc' ):
        nt,nzE,ncol = np.shape( a_zERA )
        nt,nz, ncol = np.shape( a_zCAM )
        for i in np.arange( nt ):
            for c in np.arange( ncol ):
                zoo=np.where( pmid_zCAM[i,:,c] > ps_ERA[i,c] )
                lzoo=len( zoo )
                for z in np.arange( start=0, stop=lzoo, step=1 ):
                    a_zCAMf[i,zoo[z],c] = a_zERA[i,nzE-1,c]
                    
    if ( Gridkey == 'tzyx' ):
        nt,nzE,ny,nx = np.shape( a_zERA)
        nt,nz, ny,nx = np.shape( a_zCAM )
        for i in np.arange( nt ):
            for y in np.arange( ny ):
                for x in np.arange( nx ):
                    zoo=np.where( pmid_zCAM[i,:,y,x] > ps_ERA[i,y,x] )
                    lzoo=len( zoo )
                    for z in np.arange( start=0, stop=lzoo, step=1 ):
                        a_zCAMf[i,zoo[z],y,x] = a_zERA[i,nzE-1,y,x]
                        
    toc = time.perf_counter()
    ProcTime = f" ... Bottom filling took {toc - tic:0.4f} seconds"
    print(ProcTime)
        

    return a_zCAMf



def PsAdjust( phis, phis_CAM, ps, pm150, te150, Gridkey ):

    print( " In PsAdjust " )
    print( 'PHIS', np.shape( phis ) )
    print( 'PHIS_CAM', np.shape( phis_CAM ) )
    print( 'PS', np.shape( ps ) )

    dTdz = -6.5 * 1.e-3
    t_ref1    = 290.5
    t_ref2    = 255.0
    threshold = 0.001


    lapse = -dTdz

    # Shape independent calcs
    #---------------------------
    del_phis = phis - phis_CAM # shape = H
    tsurf = te150 * ( 1.0 + lapse * (Rdry/grav) * ( ( ps/pm150 ) - 1. ) ) # shape = TH
    
    if ( Gridkey == 'tzyx' ):
        nt,ny,nx = np.shape( ps )
        ps_new = np.zeros(  (nt,ny,nx) )
        for i in np.arange( nt ):
            for y in np.arange( ny ):
                for x in np.arange( nx ):
                    t0    = tsurf[i,y,x] + lapse*phis[y,x]/grav
                    
                    #if (t0 .gt. t_ref1 .and. tsurf[i,y,x] .le. t_ref1):
                    if ( (t0 > t_ref1) and (tsurf[i,y,x] <= t_ref1 ) ):
                        lapse = ( t_ref1 - tsurf[i,y,x] )*grav/phis[y,x]

                    #elif (t0 .gt. t_ref1 .and. tsurf[i,y,x] .gt. t_ref1):
                    elif ( (t0   >   t_ref1)  and  (tsurf[i,y,x]  >   t_ref1) ):
                        lapse = 0.
                        tsurf[i,y,x] = (t_ref1 + tsurf[i,y,x] )*0.5
                    

                    #if(tsurf[i,y,x] .lt. t_ref2):
                    if(tsurf[i,y,x]   <   t_ref2):
                        lapse = -dTdz
                        tsurf[i,y,x] = (t_ref2 + tsurf[i,y,x] )*0.5


                    xx   = lapse* del_phis[y,x] / ( grav*tsurf[i,y,x]  )
                    tmp = 1.0 - xx/2.0 + xx**2.0/3.0
                    tmp = ( del_phis[y,x]  /(  Rdry*tsurf[i,y,x] ) ) *tmp
                    tmp_simple = ( del_phis[y,x]  /(  Rdry*te150[i,y,x] ) )
                    ps_new[i,y,x] = ( ps[i,y,x] ) *np.exp( tmp_simple )


    if ( Gridkey == 'tzc' ):
        nt,ncol = np.shape( ps )
        ps_new = np.zeros(  (nt,ncol ) )
        for i in np.arange( nt ):
            for c in np.arange( ncol ):
                t0    = tsurf[i,c] + lapse*phis[c]/grav
                    
                #if (t0 .gt. t_ref1 .and. tsurf[i,c] .le. t_ref1):
                if (t0    >  t_ref1  and  tsurf[i,c]  <=  t_ref1):
                    lapse = (t_ref1 - tsurf[i,c] )*grav/phis[c]

                #elif (t0 .gt. t_ref1 .and. tsurf[i,c] .gt. t_ref1):
                elif (t0   >   t_ref1  and  tsurf[i,c]   >  t_ref1):
                    lapse = 0.
                    tsurf[i,c] = (t_ref1 + tsurf[i,c] )*0.5
                    

                #if(tsurf[i,c] .lt. t_ref2):
                if(tsurf[i,c]   <   t_ref2):
                    lapse = -dTdz
                    tsurf[i,c] = (t_ref2 + tsurf[i,c] )*0.5


                xx   = lapse*del_phis[c]/(grav*tsurf[i,c]  )
                tmp = 1. - xx/2. + xx**2./3.
                tmp = del_phis[c]/(  Rdry*tsurf[i,c] )*tmp
                ps_new[i,c] = ps[i,c]*np.exp(tmp)


    return ps_new
