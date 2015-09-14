
""" A python version of the specsynth FORTRAN code.
Works entirely as a funciton
"""


import numpy as np
from scipy import signal as sg
import files as fi
import model as mo



def spectrum(zabs, zcell, logN, dopplerb, cellID, ion, vmax, inst, ntrans, 
              transName, lamb0, fosc, gamma):

    """
    Funtional form of specsynth.f

    Accepts:
        zabs (redshift of absorption)
        zcell (array of cell redshifts)
        logN (array of cell column densities)
        dopplerb (array of cell Doppler b parameters)
        cellID (array of cell ID's)
        ion (name of ion, such as CIV or MgII)
        vmax (maximum velocity of absorption window)
        inst (name of instrument, such as COSNUV or HIRES)
        transName (name of transition, such as CIV1548)
        lamb0 (central wavelength of absorption)
        fosc (oscillator strength of transition)
        gamma (damping constant of transition)

    Returns:
        lamb (array of wavelengths of spectrum)
        flux (array of fluxes of spectrum)
        
    """

    # Open line list file
#        readlines( linelist[i], trani[i] )
    # Nope, this is passed in

    # Assign atomic constants to all lines for this transition and 
    # obtain this transitions intrumetnal configuration
    con1, con2 = mo.set_atomic( lamb0, fosc, gamma )
    
    spec =  fi.config_spec( instrument, vmax, zabs, lamb0) #instrlist, j, wcen )
    R_isf, presel, rn, waveMin, waveMax, wcen, dwave = spec

    # Initialize the spectrum continuum
    ndata, lamb, wrkflux = mo.init_sectrum(waveMin, waveMax, dwave)

    m = ndata

    # Setup the ISF FSWM over the spectrum
    mo.instrument(m, wcen, 0)
    
    # Line flux radiative transfer
    wrkflux = do_abs_lines(lamb0, zabs, nline, zline, bline, 
                           con1, con2, lamb, wrkflux)

    # Lyman limit break
    dolymanlimit

    # Convolve with ISF 
    if (convolving):
        colvolve

    
    # Return the spectrum
    return velocity, flux










            






