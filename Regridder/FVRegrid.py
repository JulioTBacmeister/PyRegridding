#import xarray as xr
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGi

def uvStaggers( U,V,lat,lon ):
    nt,nz,ny,nx=np.shape(U)
    lonx = np.concatenate( ( lon[nx-1].reshape(1) ,lon), axis=0)
    lonx[0]=lonx[0]-360.
    
    slon=(lonx[1:]+lonx[0:-1] )/2.
    slat=(lat[1:]+lat[0:-1] )/2.
    VS=np.zeros( (nt, nz, ny, nx) )
    US=np.zeros( (nt, nz, ny-1, nx) )
    
    VLat, VLon = np.meshgrid(lat, slon )
    for t in np.arange( nt ):
        for z in np.arange( nz ):
            V_l   = V[ t, z, :, : ]
            V_lx  = np.concatenate( ( V_l[:,nx-1].reshape(ny,1), V_l) , axis=1)
            intrV = RGi( (lat,lonx),  V_lx )
            VS_l  = intrV( (VLat, VLon) ).T
            VS[t,z,:,:] = VS_l
            
    ULat, ULon = np.meshgrid(slat, lon )
    for t in np.arange( nt ):
        for z in np.arange( nz ):
            U_l   = U[ t, z, :, : ]
            intrV = RGi( (lat,lon),  U_l )
            US_l  = intrV( (ULat, ULon) ).T
            US[t,z,:,:] = US_l
                
    return US,VS,slat,slon