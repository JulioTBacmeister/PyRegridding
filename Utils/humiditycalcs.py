import xarray as xr
import numpy as np

epsilon = 0.622 # molecular wieght ratio (M_wv/M_dry) 

#------------------calculate vapor pressure from specific humidity and surface pressure
#https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
def calcvpfromhuss(huss, ps):
    """Calculate vapor pressure (in hPa) from specific humidity (in kg/kg) and 
       surface pressure (in Pa) i
    https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
    """
    e = (huss*ps)/(epsilon + 0.378*huss)
    e = e/100.
    e = e.rename('vp')
    return e
#---------------------------------------------------------------------------------------

#------------------calculate saturation vapor pressure from temperature 
#                  (or vapor pressure from dew point temperature)
#  Based on Bolton (1980) The computation of equivalent potential temperature, MWR
#-------------------------------------------------------------------------------------
def esat(TK):
    """calculate the saturation vapor pressure (in Pa) from T (in K)"""
    T = TK-273.15
    e = 611.2*np.exp( (17.67*T)/(T + 243.5))
    return e

def qsat(p,T):
    # saturation specifi humidity from T(K) and P(Pa):
    e=esat(T)
    qs = epsilon * e /( p - (1-epsilon)* e )
    return qs
