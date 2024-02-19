import os
#-----------------------------------------
# Find path to this module 
#------------------------------------------
module_dir = os.path.dirname(os.path.abspath(__file__))
print( f" {__name__} in {module_dir} ")



def pi():
    value=3.141592653589793238
    units='1'
    return value
def Rdry():
    value= 287.05
    units='J K-1 kg-1 '
    return value
def Rdair():
    value= avogadro()*boltzmann()/mwdair()
    units='J K-1 kg-1 '
    return value
def grav():
    value= 9.80616
    units='m s-2'
    return value
def boltzmann():
    value = 1.38065e-23
    units = 'J K-1 molecule-1'
    return value
def avogadro():
    value = 6.02214e26
    units = "molecule kmole-1"
    return value
def mwdair():
    value = 28.966
    units = "kg kmole-1"
    return value
def cpair():
    value= (7./2.)* avogadro()*boltzmann()/mwdair() # 1005.0
    units='J K-1 kg-1 '
    return value
def cpdair():
    # 3 translation + 2 rotational dof + 
    value= (7./2.)* avogadro()*boltzmann()/mwdair() # 1005.0
    units='J K-1 kg-1 '
    return value
def cvdair():
    # 3 translation + 2 rotational dof
    value= (5./2.)* avogadro()*boltzmann()/mwdair() 
    units='J K-1 kg-1 '
    return value
def Rearth():
    # Earth Radius in meters
    value= 6_378_137.0
    units='m'
    return value



"""
    gravit  = 9.80616           # acceleration of gravity ~ m/s^2
    boltz   = 1.38065e-23       # boltzmann's constant ~ J/k/molecule
    avogad  = 6.02214e26        # avogadro's number ~ molecules/kmole
    mwdair  = 28.966            # molecular weight dry air ~ kg/kmole

    rgas    = avogad*boltz      # universal gas constant ~ J/k/kmole
    rdair   = rgas/mwdair       #constant for dry air   ~ J/k/kg
"""