
""" A python version of the specsynth FORTRAN code.
Works entirely as a funciton
"""


import numpy as np
from scipy import signal as sg
import files as fi
import model as mo



def spectrum(zabs, zline, nline, bline, cellID, ion, vmax, inst, 
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
    
    spec =  fi.config_spec( inst, vmax, zabs, lamb0) #instrlist, j, wcen )
    R_isf, presel, rn, waveMin, waveMax, wcen, dwave = spec

    # Initialize the spectrum continuum
    ndata, lamb, wrkflux = mo.init_spectrum(waveMin, waveMax, dwave)

    m = ndata

    # Setup the ISF FSWM over the spectrum
    response =  mo.instrument(m, wcen, 0, R_isf, dwave)

    # Line flux radiative transfer
    wrkflux = mo.do_abs_lines(lamb0, zabs, nline, zline, bline, 
                           con1, con2, lamb, wrkflux)


    # Convolve with ISF 
    rawflux = wrkflux
    flux = sg.convolve(wrkflux, response)

    # Get the velocity of the spectrum
    velocity = mo.calc_velocity(lamb, zabs, lamb0)

    # print to file
    f = open('test.spec', 'w')
    for i in range(0,len(lamb)):
        f.write('{0:.6f}\t{1:.6f}\t{2:.6f}\n'.format(lamb[i], velocity[i], rawflux[i]))
    f.close()
    f = open('test.convolve', 'w')
    for i in range(0,len(flux)):
        f.write('{0:f}\n'.format(flux[i]))
    f.close()
    f = open('isf.dat', 'w')
    for i in range(0,len(response)):
        f.write('{0:f}\n'.format(response[i]))
    f.close()


    # Return the spectrum
    return lamb, velocity, wrkflux










            






