#!/usr/bin/env python
# Import packages 
import sys
import os 

import xarray as xr
import numpy as np
import pandas as pd
import scipy

import FVStagger as FV

import time as TimeUtil

#from GlobalVarClass import Gv as *
from GlobalVarClass import Gv



def write_netcdf( version='' ):

    
    pdTime_ERA = Gv.pdTime_ERA
    
    aint_CAM = Gv.aint_CAM
    bint_CAM = Gv.bint_CAM    
    amid_CAM = Gv.amid_CAM
    bmid_CAM = Gv.bmid_CAM
    
    lon_CAM = Gv.lon_CAM    
    lat_CAM = Gv.lat_CAM
    
    area_CAM = Gv.area_CAM
    phis_CAM = Gv.phis_CAM
    phis_ERA_xCAM = Gv.phis_ERA_xCAM
    ps_CAM = Gv.ps_CAM
    
    te_ERA_xzCAM = Gv.te_ERA_xzCAM
    q_ERA_xzCAM  = Gv.q_ERA_xzCAM
    u_ERA_xzCAM  = Gv.u_ERA_xzCAM
    v_ERA_xzCAM  = Gv.v_ERA_xzCAM
    w_ERA_xzCAM  = Gv.w_ERA_xzCAM
    
    if (Gv.MySrc=="ERA5"):
        SuperDir = "/glade/campaign/cgd/amp/juliob/ERA5"
    elif (Gv.MySrc=="ERAI"):
        SuperDir = "/glade/campaign/cgd/amp/juliob/ERAI"
    else:
        SuperDir = "/glade/campaign/cgd/amp/juliob/MiscRegridding"

    ntime = np.shape(pdTime_ERA)[0]
    print(ntime)
    #Bdiro="/glade/derecho/scratch/juliob/ERA5/" + Gv.MyDst    
    Bdiro=f"{SuperDir}/{Gv.MyDst}/{Gv.MyDstVgrid}"
    #######
    os.makedirs( Bdiro , exist_ok=True )

    #version='test_netcdf4_default'
    #Bfilo="/glade/derecho/scratch/juliob/ERA5/" + Gv.MyDst + "/" + Gv.MySrc +"_x_"+ Gv.MyDst + "_"+ Gv.MyDstVgrid + "_" + version 
    Bfilo= Bdiro + "/" + Gv.MySrc +"_x_"+ Gv.MyDst + "_"+ Gv.MyDstVgrid + "_" + version 

    
    if (Gv.doWilliamsonOlson==True):
        Bfilo = Bfilo + '_WO'

    ilev = (aint_CAM+bint_CAM)* 1_000. #* 100_000.
    lev  = (amid_CAM+bmid_CAM)* 1_000. #* 100_000.



    if (Gv.dstTZHkey == 'tzc' ):
        nt,nz,ncol = np.shape( te_ERA_xzCAM )
        for itim in np.arange( ntime ):
            dims   = ["ncol","time","lev","ilev"]
            coords = dict( 
                lon  = ( ["ncol"],lon_CAM ),
                lat  = ( ["ncol"],lat_CAM ),
                lev  = ( ["lev"],lev),
                ilev = ( ["ilev"],ilev),
                time = ( ["time"],  np.array(itim ,ndmin=1 ,dtype=np.int32 ) ), #pd.to_datetime( pdTime_ERA[itim] ) ),
            )
        
            Wds = xr.Dataset( coords=coords  )
            #Wds["TimeStamp"] = np.array( [itim * 24./ntime] ).astype( np.int32)  # pd.to_datetime( pdTime_ERA[itim] )
            Wds["TimeStamp"] = pd.to_datetime( pdTime_ERA[itim] )
            Wds["P_00"] = 100_000.
        
            Dar = xr.DataArray( data=aint_CAM, dims=('ilev',),
                                attrs=dict( description='interface hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyai'] = Dar

            Dar = xr.DataArray( data=bint_CAM, dims=('ilev',),
                                attrs=dict( description='interface hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybi'] = Dar

            Dar = xr.DataArray( data=amid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyam'] = Dar

            Dar = xr.DataArray( data=bmid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybm'] = Dar
        
            Dar = xr.DataArray( data=area_CAM, dims=('ncol',),
                                attrs=dict( description='Cell area',units='Steradians',) ,) 
            Wds['area'] = Dar

            Dar = xr.DataArray( data=phis_CAM, dims=('ncol',),
                                attrs=dict( description='Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS'] = Dar

            Dar = xr.DataArray( data=phis_ERA_xCAM, dims=('ncol',),
                                attrs=dict( description='ERA Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS_ERA'] = Dar

            Dar = xr.DataArray( data=ps_CAM[itim,:].reshape(1,ncol), 
                                dims=('time','ncol',),
                                attrs=dict( description='Surface Pressure',units='Pa',) ,) 
            Wds['PS'] = Dar
    
            Dar = xr.DataArray( data=te_ERA_xzCAM[itim,:,:].reshape(1,nz,ncol), 
                                dims=('time','lev','ncol',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['T'] = Dar

            Dar = xr.DataArray( data=q_ERA_xzCAM[itim,:,:].reshape(1,nz,ncol), 
                                dims=('time','lev','ncol',),
                                attrs=dict( description='specific humidity',units='kg kg-1',) ,) 
            Wds['Q'] = Dar
        
            Dar = xr.DataArray( data=u_ERA_xzCAM[itim,:,:].reshape(1,nz,ncol), 
                                dims=('time','lev','ncol',),
                                attrs=dict( description='X-wind',units='m s-1',) ,) 
            Wds['U'] = Dar

            Dar = xr.DataArray( data=v_ERA_xzCAM[itim,:,:].reshape(1,nz,ncol), 
                                dims=('time','lev','ncol',),
                                attrs=dict( description='Y-wind',units='m s-1',) ,) 
            Wds['V'] = Dar

            Dar = xr.DataArray( data=w_ERA_xzCAM[itim,:,:].reshape(1,nz,ncol), 
                                dims=('time','lev','ncol',),
                                attrs=dict( description='Vertical motion',units='Pa s-1',) ,) 
            Wds['W'] = Dar

            yymmdd = str(pdTime_ERA[itim])[0:10]
            hr=str(pdTime_ERA[itim])[11:13]
            ss = str(int(hr)*3600).zfill(5)
            timetag =  yymmdd+'-'+ss
            filo= Bfilo + "." + timetag+ ".nc"
            print( filo )
            Wds.to_netcdf( filo  ) #, format="NETCDF3_64BIT" , engine='scipy' ) #,format="NETCDF3_CLASSIC" )

    if (Gv.dstTZHkey == 'tzyx' ):
        tic_FVstag = TimeUtil.perf_counter()
        US,VS,slat,slon = FV.uvStaggers(U=u_ERA_xzCAM, 
                                        V=v_ERA_xzCAM,
                                        lon=lon_CAM,
                                        lat=lat_CAM   )
        toc_FVstag = TimeUtil.perf_counter()
        pTime = f"Creating FV staggered US,VS took {toc_FVstag - tic_FVstag:0.4f} seconds"
        print(pTime)
        
        nt,nz,ny,nx = np.shape( te_ERA_xzCAM )
        for itim in np.arange( ntime ):
            dims   = ["lon","lat","time","lev","ilev"]
            coords = dict( 
                lon  = ( ["lon"],lon_CAM ),
                lat  = ( ["lat"],lat_CAM ),
                slon  = ( ["slon"],slon ),
                slat  = ( ["slat"],slat ),
                lev  = ( ["lev"],lev),
                ilev = ( ["ilev"],ilev),
                time = ( ["time"],  np.array(itim ,ndmin=1 ) ), #pd.to_datetime( pdTime_ERA[itim] ) ),
            )
        
            Wds = xr.Dataset( coords=coords  )
            Wds["TimeStamp"] = pd.to_datetime( pdTime_ERA[itim] )
            Wds["P_00"] = 100_000.
        
            Dar = xr.DataArray( data=aint_CAM, dims=('ilev',),
                                attrs=dict( description='interface hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyai'] = Dar

            Dar = xr.DataArray( data=bint_CAM, dims=('ilev',),
                                attrs=dict( description='interface hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybi'] = Dar

            Dar = xr.DataArray( data=amid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate A-coeff ',units='1',) ,) 
            Wds['hyam'] = Dar

            Dar = xr.DataArray( data=bmid_CAM, dims=('lev',),
                                attrs=dict( description='mid-level hybrid eta coordinate B-coeff ',units='1',) ,) 
            Wds['hybm'] = Dar
        
            Dar = xr.DataArray( data=area_CAM, dims=('lat','lon',),
                                attrs=dict( description='Cell area',units='Steradians',) ,) 
            Wds['area'] = Dar

            Dar = xr.DataArray( data=phis_CAM, dims=('lat','lon',),
                                attrs=dict( description='Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS'] = Dar

            Dar = xr.DataArray( data=phis_ERA_xCAM, dims=('lat','lon',),
                                attrs=dict( description='ERA Surface Geopotential Height',units='m+2 s-2',) ,) 
            Wds['PHIS_ERA'] = Dar

            Dar = xr.DataArray( data=ps_CAM[itim,:,:].reshape(1,ny,nx) , 
                                dims=('time','lat','lon',),
                                attrs=dict( description='Surface Pressure',units='Pa',) ,) 
            Wds['PS'] = Dar
    
            Dar = xr.DataArray( data=te_ERA_xzCAM[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','lon',),
                                attrs=dict( description='Air Temperature',units='K',) ,) 
            Wds['T'] = Dar

            Dar = xr.DataArray( data=q_ERA_xzCAM[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','lon',),
                                attrs=dict( description='Specific humidity',units='kg kg-1',) ,) 
            Wds['Q'] = Dar
        
            Dar = xr.DataArray( data=u_ERA_xzCAM[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','lon',),
                                attrs=dict( description='X-wind',units='m s-1',) ,) 
            Wds['U'] = Dar

            Dar = xr.DataArray( data=v_ERA_xzCAM[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','lon',),
                                attrs=dict( description='Y-wind',units='m s-1',) ,) 
            Wds['V'] = Dar

            Dar = xr.DataArray( data=US[itim,:,:,:].reshape(1,nz,ny-1,nx), 
                                dims=('time','lev','slat','lon',),
                                attrs=dict( description='Staggered X-wind',units='m s-1',) ,) 
            Wds['US'] = Dar

            Dar = xr.DataArray( data=VS[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','slon',),
                                attrs=dict( description='Staggered Y-wind',units='m s-1',) ,) 
            Wds['VS'] = Dar

            Dar = xr.DataArray( data=w_ERA_xzCAM[itim,:,:,:].reshape(1,nz,ny,nx), 
                                dims=('time','lev','lat','lon',),
                                attrs=dict( description='Vertical motion',units='Pa s-1',) ,) 
            Wds['W'] = Dar
        
            yymmdd = str(pdTime_ERA[itim])[0:10]
            hr=str(pdTime_ERA[itim])[11:13]
            ss = str(int(hr)*3600).zfill(5)
            timetag =  yymmdd+'-'+ss
            filo= Bfilo + "." + timetag+ ".nc"
            print( filo )
            Wds.to_netcdf( filo ) #,format="NETCDF3_CLASSIC" )

    code = 1
    return code

def ExFromGv():
    
    pdTime_ERA = Gv.pdTime_ERA
    
    
    aint_CAM = Gv.aint_CAM
    bint_CAM = Gv.bint_CAM    
    amid_CAM = Gv.amid_CAM
    bmid_CAM = Gv.bmid_CAM
    
    lon_CAM = Gv.lon_CAM    
    lat_CAM = Gv.lat_CAM
    
    area_CAM = Gv.area_CAM
    phis_CAM = Gv.phis_CAM
    phis_ERA_xCAM = Gv.phis_ERA_xCAM
    ps_CAM = Gv.ps_CAM
    
    te_ERA_xzCAM = Gv.te_ERA_xzCAM
    q_ERA_xzCAM  = Gv.q_ERA_xzCAM
    u_ERA_xzCAM  = Gv.u_ERA_xzCAM
    v_ERA_xzCAM  = Gv.v_ERA_xzCAM
    w_ERA_xzCAM  = Gv.w_ERA_xzCAM
    
    
    
    
    
