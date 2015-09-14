# Functions for creating the model absorption for spectrum.py

import numpy as np
from scipy import constants as sc
from astropy import constants as const

def setatomic(j, lam0, f0):

    c = const.c.to('cm/s')

    con1Const = ( np.sqrt(sc.pi)*np.pow(e,2) ) / ( sc.m_e*np.pow(sc.c, 2) )
    con2Const = 4 * sc.pi * sc.c
    con1 = 1e5*f0*np.pow(lam0,2) * con1Const


    return con1, con2

