#!/usr/bin/env python
# Import packages 
import sys
import argparse as arg

import GlobalVarClass
from GlobalVarClass import Gv


import xarray as xr
import numpy as np
import pandas as pd

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


# "ChatGPI version" --- 
import VertRegridFlexLL as vrg
print( "Using Flexible parallel/serial VertRegrid ")

# import modules in other directories
sys.path.append('../Utils/')
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

#importlib.reload( Globals2 )

# Physical Constants
Rdry = Con.Rdry() # 


#Gv = Globals2.VariableContainer()

#-------------------------------------------------------------
#  Naming conventions
#-------------------------------------------------------------
# aaa_{CAM,ERA} 
# Indicates the immediate provenance of a variable, e.g.,
#      phis_CAM ==> phis from CAM on CAM's grid
#      phis_ERA ==> phis from ERA on ERA's grid
# lower case 'phis' indicates this is ndarray-like 
#
# aaa_{CAM,ERA}_x{ERA,CAM}
# Indicates variable has been remapped horizontall. So, e.g.
#      phis_ERA_xCAM ==> ERA phis remapped horizontally to the CAM grid 
#
# aaa_{CAM,ERA}_xz{ERA,CAM}
# Indicates variable has been remapped horizontally AND vertically. So, e.g.
#      te_ERA_xzCAM ==> ERA temperature remapped horizontally to the CAM horizontal grid 
#                       and then also vertically interpoated to the CAM vertical grid
#-------------------------------------------------------------
def xRegrid( ExitAfterTemperature=False , 
             HorzInterpLnPs=False , 
             Use_ps_ERA_xCAM_in_vert=True ):


    StartTime = time.asctime( time.localtime(time.time()) )
    tic_overall = time.perf_counter()
    print( f"starting xRegrid {Gv.MySrc} _x_ {Gv.MyDst} at {StartTime} ")
    #-----------------------------------------
    # Horz remap of PHIS - No time dimension
    #-----------------------------------------
    Gv.phis_ERA_xCAM = erg.HorzRG( aSrc = Gv.phis_ERA , 
                                regrd = Gv.regrd , 
                                srcField=Gv.srcf , 
                                dstField=Gv.dstf , 
                                srcGridkey=Gv.srcHkey ,
                                dstGridkey=Gv.dstHkey )
    
    toc_here = time.perf_counter()
    pTime = f"Finished phis Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)

    """
    #------------------------------------------
    # Step in with Islas's 1440x720 =>xCAM 
    # f09 regridded ERA5 topo.
    # This is for diagnostics/debugging only
    #-----------------------------------------
    fTopo2='/glade/u/home/islas/for/suqin/regridera5/makephis/output/PHIS_model_and_ERA5_analysis_f09_f09.nc'
    dTopo2=xr.open_dataset( fTopo2 )
    phis_ERA_xCAM =dTopo2['PHIS_analysis'].values[0,:,:]
    """
 
    
    #-----------------------------------------
    # Calculate difference between phis's
    #-----------------------------------------
    Dphis = Gv.phis_ERA_xCAM - Gv.phis_CAM
    
        
    #-----------------------------------------
    # Horz remap of PS
    # log-exp bit is to reproduce W&O
    #-----------------------------------------
    if (HorzInterpLnPs==True):
        xfld_ERA = np.log( Gv.ps_ERA )
    else:
        xfld_ERA = Gv.ps_ERA
        
    Gv.ps_ERA_xCAM    = erg.HorzRG( aSrc = xfld_ERA , 
                                 regrd = Gv.regrd , 
                                 srcField=Gv.srcf , 
                                 dstField=Gv.dstf , 
                                 srcGridkey=Gv.srcTHkey ,
                                 dstGridkey=Gv.dstHkey )

    if (HorzInterpLnPs==True):
        Gv.ps_ERA_xCAM = np.exp( Gv.ps_ERA_xCAM ) 

    toc_here = time.perf_counter()
    pTime = f"Finished ps Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)

    #-----------------------------------------
    # Make 3D prssure field on ERA ZH grid
    #-----------------------------------------
    Gv.pmid_ERA, Gv.pint_ERA, Gv.delp_ERA \
        = MkP.Pressure (am=Gv.amid_ERA ,
                        bm=Gv.bmid_ERA ,
                        ai=Gv.aint_ERA ,
                        bi=Gv.bint_ERA ,
                        ps=Gv.ps_ERA ,
                        p_00=Gv.p_00_ERA ,
                        Gridkey = Gv.srcTZHkey )

    #-----------------------------------------
    # Horz remap of temperature
    #-----------------------------------------
    Gv.te_ERA_xCAM    = erg.HorzRG( aSrc = Gv.te_ERA , 
                                 regrd = Gv.regrd , 
                                 srcField=Gv.srcf , 
                                 dstField=Gv.dstf , 
                                 srcGridkey=Gv.srcTZHkey ,
                                 dstGridkey= Gv.dstHkey )
    
    # Save off a copy of te before any funny business 
    # happens, e.g., Williamson&Olson
    #------------------------------------------------
    Gv.te_ERA_xCAM_00 = copy.deepcopy( Gv.te_ERA_xCAM  )
    toc_here = time.perf_counter()
    pTime = f"Finished te_ERA Horz Rgrd  {toc_here - tic_overall:0.4f} seconds"
    print(pTime)

    #-----------------------------------------------------------
    # Find "T_bot" and "P_bot" defined in Williamson & Olson 
    # as occurring at "the first level above 150m ... " 
    # This uses ERA temperature and surface pressure remapped  
    # to CAM horz grid (adapted from Isla's function)
    #-----------------------------------------------------------
    
    tic_150 = time.perf_counter()

    Gv.te_150, Gv.pmid_150, Gv.L150  =\
                    MkP.Pressure_TandP150(
                    am=Gv.amid_ERA ,
                    bm=Gv.bmid_ERA ,
                    ai=Gv.aint_ERA ,
                    bi=Gv.bint_ERA ,
                    ps=Gv.ps_ERA_xCAM ,
                    te=Gv.te_ERA_xCAM , 
                    p_00=Gv.p_00_ERA , 
                    Gridkey = Gv.dstTZHkey )


    toc_150 = time.perf_counter()
    pTime = f"Finding Te150 and P150 took  {toc_150 - tic_150:0.4f} seconds"
    print(pTime)

    
    #-------------------------------------------------------------------
    #                    "CAM surface pressure"
    #-------------------------------------------------------------------
    # We don't actually have ps from CAM, so we make a guess based on
    # ps_ERA_xCAM and te_bot asdescribed above.  In a sense this is a 
    # vertical remapping, so we could call this variable ps_ERA_xzCAM, 
    # but we'll just call it ps_CAM ...
    #-------------------------------------------------------------------
    
    
    Gv.ps_CAM = vrg.PsAdjust( phis=Gv.phis_ERA_xCAM, 
                            phis_CAM=Gv.phis_CAM, 
                            ps= Gv.ps_ERA_xCAM , 
                            pm150=Gv.pmid_150 , 
                            te150=Gv.te_150 , 
                            Gridkey=Gv.dstTZHkey  )

   
    #-----------------------------------------------------------------------------------------------
    # Now we creat full 4(3)D pressure fields on the ERA and CAM vertical grids. These are used for 
    # vertical interpolation below. Not clear what surface pressure to use when building ERA vertical
    # grid, i.e., ps_CAM or ps_ERA_xCAM.
    #
    # The WO2015 document seems to suggest ps_CAM but W&O code definitely uses ps_ERA_xCAM
    #-----------------------------------------------------------------------------------------------
    tic_P3D = time.perf_counter()
    
    # Surface pressure: Choose wisely
    #---------------------------------
    if (Use_ps_ERA_xCAM_in_vert == True):
        ps_FLD = Gv.ps_ERA_xCAM
    else:
        ps_FLD = Gv.ps_CAM
        
    Gv.pmid_CAM_zERA, Gv.pint_CAM_zERA, Gv.delp_CAM_zERA \
        = MkP.Pressure (am=Gv.amid_ERA ,
                        bm=Gv.bmid_ERA ,
                        ai=Gv.aint_ERA ,
                        bi=Gv.bint_ERA ,
                        ps=ps_FLD , # What to use here seems key: ps_CAM or ps_ERA_xCAM
                        p_00=Gv.p_00_ERA ,
                        Gridkey = Gv.dstTZHkey )

    Gv.pmid_CAM,Gv.pint_CAM,Gv.delp_CAM \
        = MkP.Pressure (am=Gv.amid_CAM ,
                        bm=Gv.bmid_CAM ,
                        ai=Gv.aint_CAM ,
                        bi=Gv.bint_CAM ,
                        ps=Gv.ps_CAM ,
                        p_00=Gv.p_00_CAM , 
                        Gridkey = Gv.dstTZHkey )

    #-----------------------------------------------------------
    # Log-pressure is preferred for vertical interpolation
    # per Williamson&Olson
    #-----------------------------------------------------------
    p_00 = 100_000. # Here we just use the sensible value of p_00
    lnpint_CAM = -7_000. * np.log( Gv.pint_CAM / p_00 )
    lnpmid_CAM = -7_000. * np.log( Gv.pmid_CAM / p_00 )
    lnpmid_CAM_zERA = -7_000. * np.log( Gv.pmid_CAM_zERA /p_00 )

    toc_P3D = time.perf_counter()
    pTime = f"Creating 3D P-fields etc., took   {toc_P3D - tic_P3D:0.4f} seconds"
    print(pTime)
    
    
    if ( Gv.doWilliamsonOlson == True ):
        tic_WO = time.perf_counter()
        print( "WilliamsonOlson surface " )
        #----------------------------------------------------------
        # Calculate extrapolated surface temperature using 
        # Williamson & Olson standard lapse rate approach
        #----------------------------------------------------------
        Gv.ts_extrap = vrg.TsExtrap( ps = Gv.ps_CAM ,
                                  pm150 = Gv.pmid_150 ,
                                  te150 = Gv.te_150 )

                              
        #--------------------------------------------------------
        # If Williamson & Olson treatment of surface layer 
        # is selected then correct te_ERA_xzCAM between 
        # pmid_150 and ps_CAM
        #-------------------------------------------------------
        Gv.te_WO  =    vrg.TeWO( te = Gv.te_ERA_xCAM ,
                              pmid = Gv.pmid_CAM_zERA ,
                              te150 = Gv.te_150 ,
                              pm150 = Gv.pmid_150 ,
                              ts = Gv.ts_extrap,
                              ps = Gv.ps_CAM, 
                              L150 = Gv.L150 ,
                              Gridkey = Gv.dstTZHkey )

        Gv.te_ERA_xCAM = copy.deepcopy( Gv.te_WO )
        toc_WO = time.perf_counter()
        pTime = f"Williamson Olson surface took  {toc_WO - tic_WO:0.4f} seconds"
        print(pTime)                       
    

    print(" going into vertical regrid of T " )
    Gv.te_ERA_xzCAM = vrg.VertRG( a_x  = Gv.te_ERA_xCAM ,
                               zSrc = lnpmid_CAM_zERA ,
                               zDst = lnpmid_CAM ,
                               Gridkey =Gv.dstTZHkey ,
                               kind = 'quadratic' ) #linea
    

    #-------------------------------------------------------------
    # return statement to assist in debugging and analysis
    #-------------------------------------------------------------
    if ( ExitAfterTemperature == True ):
        return


    #-------------------------------------------------------------
    # Continue on to regridding Q, U, V, W ...
    #-------------------------------------------------------------

    #--------------------
    #  Regridding of Q
    #---------------------
    print(" going into horz+vertical regrid of Q " )
    Gv.q_ERA_xzCAM , Gv.q_ERA_xCAM =\
                        fullRegrid( a_ERA = Gv.q_ERA ,
                                    zSrc = lnpmid_CAM_zERA ,
                                    zDst = lnpmid_CAM )
        
    Gv.q_ERA_xzCAM = vrg.BottomFill( a_zCAM = Gv.q_ERA_xzCAM ,
                                  a_zERA = Gv.q_ERA_xCAM ,
                                  pmid_zCAM=Gv.pmid_CAM ,
                                  ps_ERA = Gv.ps_ERA_xCAM , 
                                  Gridkey = Gv.dstTZHkey )

    qx = SaturateQ( q=Gv.q_ERA_xzCAM , 
                    te=Gv.te_ERA_xzCAM ,
                    p=Gv.pmid_CAM, 
                    Gridkey = Gv.dstTZHkey )
    
    Gv.q_ERA_xzCAM =  copy.deepcopy(qx)




    #--------------------
    #  Regridding of U
    #---------------------
    print(" going into horz+vertical regrid of U " )
    Gv.u_ERA_xzCAM, Gv.u_ERA_xCAM = fullRegrid ( a_ERA = Gv.u_ERA ,
                                           zSrc = lnpmid_CAM_zERA ,
                                           zDst = lnpmid_CAM )

    Gv.u_ERA_xzCAM = vrg.BottomFill( a_zCAM = Gv.u_ERA_xzCAM ,
                                  a_zERA = Gv.u_ERA_xCAM ,
                                  pmid_zCAM=Gv.pmid_CAM ,
                                  ps_ERA = Gv.ps_ERA_xCAM , 
                                  Gridkey = Gv.dstTZHkey )
    
    #--------------------
    #  Regridding of V
    #---------------------
    print(" going into horz+vertical regrid of V " )
    Gv.v_ERA_xzCAM, Gv.v_ERA_xCAM = fullRegrid ( a_ERA = Gv.v_ERA ,
                                           zSrc = lnpmid_CAM_zERA ,
                                           zDst = lnpmid_CAM )

    Gv.v_ERA_xzCAM = vrg.BottomFill( a_zCAM = Gv.v_ERA_xzCAM ,
                                  a_zERA = Gv.v_ERA_xCAM ,
                                  pmid_zCAM=Gv.pmid_CAM ,
                                  ps_ERA = Gv.ps_ERA_xCAM , 
                                  Gridkey = Gv.dstTZHkey )
    
    #--------------------
    #  Regridding of W
    #---------------------
    print(" going into horz+vertical regrid of W " )
    Gv.w_ERA_xzCAM, Gv.w_ERA_xCAM = fullRegrid ( a_ERA = Gv.w_ERA ,
                                          zSrc = lnpmid_CAM_zERA ,
                                          zDst = lnpmid_CAM )

    Gv.w_ERA_xzCAM = vrg.BottomFill( a_zCAM = Gv.w_ERA_xzCAM ,
                                  a_zERA = Gv.w_ERA_xCAM ,
                                  pmid_zCAM=Gv.pmid_CAM ,
                                  ps_ERA = Gv.ps_ERA_xCAM , 
                                  Gridkey = Gv.dstTZHkey )
    

    
    toc = time.perf_counter()
    pTime = f"Overall time in this function  {toc - tic_overall:0.4f} seconds"
    print(pTime)
        
    rcode =1 
    return rcode

def fullRegrid( a_ERA,  zSrc ,  zDst , kind='linear', ReturnVars=2 ):
    
    print("Horz RG in fullRegrid " )
    a_ERA_xCAM    = erg.HorzRG( aSrc = a_ERA , 
                                 regrd = Gv.regrd , 
                                 srcField=Gv.srcf , 
                                 dstField=Gv.dstf , 
                                 srcGridkey=Gv.srcTZHkey,
                                 dstGridkey=Gv.dstHkey )

    print("Vert RG in fullRegrid " )
    a_ERA_xzCAM = vrg.VertRG( a_x  = a_ERA_xCAM ,
                              zSrc = zSrc ,
                              zDst = zDst ,
                              Gridkey=Gv.dstTZHkey,
                              kind = kind )

    if (ReturnVars==1):
        return a_ERA_xzCAM
    if (ReturnVars==2):
        return a_ERA_xzCAM,a_ERA_xCAM

def SaturateQ ( q , te , p, Gridkey ):

    qsat = hum.qsat( p=p, T=te )
    qx=np.minimum( q , qsat )

    """
    print( "UGLY temporary Kluge: Don't saturate/ZERO top 3 levels of model " )
    if (Gridkey == 'tzyx' ):
        qx[:,0:2,:,:]= 0*q[:,0:2,:,:]
    if (Gridkey == 'tzc' ):
        qx[:,0:2,:]= 0*q[:,0:2,:]
    if (Gridkey == 'zyx' ):
        qx[0:2,:,:]= 0*q[0:2,:,:]
    if (Gridkey == 'zc' ):
        qx[0:2,:]= 0*q[0:2,:]
    
    """
    
    return qx

