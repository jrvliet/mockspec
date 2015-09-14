# Functions for creating the model absorption for spectrum.py

import numpy as np
#from scipy import constants as sc
from astropy import constants as const

def set_atomic(lamb0, f0, gamma0):

    c = const.c.cgs.value
    e = const.e.cgs.value
    m_e = const.m_e.cgs.value    

    con1Const = ( np.sqrt(np.pi)*e*e ) / ( m_e*c*c )
    con2Const = 4 * np.pi * c 
    con1 = 1e5 * f0 * lamb0*lamb0 * con1Const
    con2 = 1e-8 * gamma0 * lamb0*lamb0 / con2Const


    return con1, con2




def init_spectrum():
    """
    Initilizes the spectrum to proper pixelization 
    and unity value
    
    Accepts:
        waveMin (minimum wavelength in the spectrum to make)
        waveMax (maximum wavelength in the spectrum to make)
        dwave (pixel sampling rate)

    Returns:
        lamb (array of wavelegth values)
        wrkflux (array of flux values, all set to 1)

    """

    # Given the wavelength range and the pixel sampling rate, 
    # determine the number of pixels "ndata" in the spectrum
    ndata = int( (waveMax - waveMin)/dwave ) + 1

    # Fill the wavelength pixel values and initialize the 
    # spectrum as unitu continuum in array wrkflux

    lamb = np.zeros(ndata)
    wrkflux = np.ones(ndata)
    for i in range(0,ndata):
        lamb[i] = waveMin + (i-1)*dwave

    return lamb, wrkflux




def do_abs_lines():


    """
    Uses Voigt profile formalism and backgorund source
    radiative transfer (no source function) to compute
    the absorption lines
    
    Some helpful definitions:
    b1  = column density
    b2  = rest frame central wavelength
    b3  = rest frame Doppler width
    y   = natural broadening in terms of doppler units
    x   = wavelength difference in doppler units
    w0  = rest frame wavelegnth for transition
    w   = rest frame wavelength fro observed wavelength
    tau = optical depth


    Accepts:
        lamb0 (central wavelength of transition)
        zabs  (redshift of absorption)
        nline (array of cell column densities)
        zline (array of cell redshifts)
        bline (array of cell Doppler b values)
        con1  (constatn 2 form set_atomic)
        con2  (constant 1 from set_atomic)
        lamb (array of wavelengths)
        wrkflux (array of flux values)

    Returns:
        wkrflux (array of flux values)

    """

    ckms = const.c.to('km/s').value

    # Initalize Voight profile parameters
    # These must be rest from quantities even though 
    # the line may be redshifted
    w0 = lamb0
    b1 = nline[i]
    b2 = w0 * (1+zline[i]) / (1+zabs)
    b3 = w0 * abs(bline[i])/ckms
    y  = con2 / b3
    tcon = con1 * (b1/b3)

    # Loop over the pixels and perfomr the radiative transfer
    for i in range(0,len(wrkflux)):
        w = lamb[i] / (1+zabs)
        x = (w-b2)/b3
        voigt(x,y,u,v)
        wrkflux[i] = wrkflux[i] * np.exp(-tcon*u)

    return wrkflux
    


def do_lyman_limit():
def drop_lines():









