
""" A python version of the specsynth FORTRAN code.
Works entirely as a funciton
"""


import numpy as np
from scipy import signal as sg



def specsynth(zabs, zcell, logN, dopplerb, cellID, ion, Vmax, inst, ntrans, 
              transName, central, fosc, gamma):

    """
    Funtional form of specsynth.f
    Takes in contents of a .lines file:
        
        zabs
        for each cell:        
            zcell       logN    Doppler b       cellID

    Also takes in the ion, Vmax, and instrument

    Returns a spectra object:
        
        Unknown         velocity        flux    unknown      unknown   unknown
        
    """

    # Open line list file
#        readlines( linelist[i], trani[i] )
    # Nope, this is passed in

    # Assign atomic constants to all lines for this transition and 
    # obtain this transitions intrumetnal configuration
    con1, con2 = setatomic( lam0, f0, gamma0 )
    
    configspec( instrlist, j, wcen )

    # Initialize the spectrum continuum
    initpsectrum 

    m = ndata

    # Setup the ISF FSWM over the spectrum
    instrument(m, wcen, 0)
    
    # Line flux radiative transfer
    doabslines

    # Lyman limit break
    dolymanlimit

    # Convolve with ISF 
    if (convolving):
        colvolve

    
    # Return the spectrum
    return velocity, flux










            






