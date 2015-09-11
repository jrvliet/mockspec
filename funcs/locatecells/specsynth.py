
""" A python version of the specsynth FORTRAN code.
Works entirely as a funciton
"""


import numpy as np
from scipy import signal as sg



def specsynth(zabs, zcell, logN, dopplerb, cellID, ion, Vmax, inst):

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

    # Get transition parameters
    paramlist = 'Mockspec.runpars'
    getparams(paramlist)
    
    # Get atomic parametes
    gettransitions(tranilist, losnum, ntran, trani, linelist)


    # Loop over the number of transistions and make the spectrum for each one
    for i in range(0,ntran):

        # Assign atomic constants to all lines for this transition and 
        # obtain this transitions intrumetnal configuration
        
    
