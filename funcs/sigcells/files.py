# Functions for manipulating files for generating
# spectra

import sys
from astropy import constants as const

def config_spec(instrument, vmax, zabs, lamb0):

    """
    Reads Mockspec.instruments and get the properties for
    the desired instrument
    
    Accepts:
        instrument (name of instrument)
        vmax (velocity range of the spectrum to make)
        zabs (redshift of the absorption feature)
        lamb0 (central wavelength of the transition)

    Returns: 
        spec (list containing:)
            R_isf (resolution)
            presel (pixels per resolution element)
            rn (read noise)
            waveMin (minimum wavelegth of absorption window)
            waveMax (maxiumum wavelength of absorption window)
            wcen (central wavelengh of absorption window)
            dwave (wavelength step size of sprectrum to make)

    """
    found = 0
    ckms = const.c.to('km/s').value

    filename = 'Mockspec.instruments'
    with open(filename) as f:
        
        # Read past header:
        f.readline()

        # Loop through file and locate the instrument name
        for line in f:
            
            instr = line.split()[0]
            if instr == instrument: 
                found = 1
                R_isf = float(line.split()[1])
                presel = float(line.split()[2])
                rn = float(line.split()[3])
                velMin = -1*vmax
                velMax = vmax

    if found == 1:

        # Compute the pixelization of the instrument
        waveMin = (1.0 + velMin/ckms) * (1.0+zabs) * lamb0
        waveMax = (1.0 + velMax/ckms) * (1.0 +zabs) * lamb0
        wcen = (waveMax + waveMin) / 2.0
        dwave = wcen / (presel*R_isf)

        spec = R_isf, presel, rn, waveMin, waveMax, wcen, dwave 
        return spec

    else:
        print 'Instrument {0:s} not found in {1:s}'.format(instrument, filename)
        print 'Exiting...'
        sys.exit()
    


    
def get_transitions(ion):
    """
    Reads in Mockspec.transitions to get the physical properties of 
    the transition the spectrum is being made for. 

    Accepts:
        ion (name of ion, such as CIV or MgII)

    Returns: 
        k (species of ion)
        j (ionization state of ion)
        transName (name of the transition)
        lamb0 (central wavelength)
        fosc (oscillator strength)
        gamma (damping constant)
        mamu (atomic mass of ion in amu)
        abund (abundance of ion)
        ip (ionization potential from ground in eV)
    """
    
    found = 0

    filename = 'Mockspec.transitions'
   
    with open(filename) as f:

        # One header line
        f.readline()

        # To be a match, column 0 (on or off) and column 4 (ion name)
        # need to match
        for line in f:
            
            l = line.split()
            
            if l[0]=='1' and l[4]==ion:
                found = 1
                k = int(l[2])
                j = int(l[3])
                transName = l[5]
                lamb0 = float(l[6])
                fosc = float(l[7])
                gamma = float(l[8])
                mamu = float(l[9])
                abund = float(l[10])
                ip = float(l[11])
                break

    if found==1:
        return k, j, transName, lamb0, fosc, gamma, mamu, abund, ip
    else:
        print 'Cannot find ion={0:s} in {1:s}'.format(ion, filename)
        sys.exit()
    







