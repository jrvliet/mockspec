# Functions for creating the model absorption for spectrum.py

import numpy as np
#from scipy import constants as sc
from astropy import constants as const

def setatomic(lam0, f0, gamma0):

    c = const.c.cgs.value
    e = const.e.cgs.value
    m_e = const.m_e.cgs.value    

    con1Const = ( np.sqrt(np.pi)*e*e ) / ( m_e*c*c )
    con2Const = 4 * np.pi * c 
    con1 = 1e5 * f0 * lam0*lamb0 * con1Const
    con2 = 1e-8 * gamma0 * lamb0*lamb0 / con2Const


    return con1, con2

