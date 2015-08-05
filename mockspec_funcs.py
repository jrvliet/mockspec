
import sys
import numpy as np


def readControlFile():

    # Reads in the control file
    # Filename: mockspec.config
    try:
        f = open('mockspec.config')
    except IOError:
        print 'Missing control file: {0:s}'.format(filename)
        print 'Quitting...'
        sys.exit()

    # Read past header
    f.readline()
    f.readline()

    galID = f.readline().split()[0]
    expn = f.readline().split()[0]
    nlos = f.readline().split()[0]
    incline = f.readline().split()[0]
    ewcut = f.readline().split()[0]
    snr = f.readline().split()[0]
    ncores = f.readline().split()[0]
    rootLoc = f.readline().split()[0]
    
    # Now at ion section
    # Loop over rest of file
    ions = []
    xh = []
    instruments = []
    f.readline()
    f.readline()
    for line in f:
        l = line.split()
        ions.append(l[0])
        xh.append(l[1])
        instruments.append(l[2])

    f.close()

    props = (galID, expn, nlos, incline, ewcut, snr, ncores, rootLoc)

    return props, ions, xh, instruments




