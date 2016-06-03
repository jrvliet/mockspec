'''
Python version of the Fortran code
'''


import pandas as pd
import numpy as np
from math import sqrt
from scipy.stats import poisson

def bindata(filename, parameter, binsize, lowerlimit, upperlimit):

    '''
    Reads in the file given by filename (assumed to be in HDF5 format)
    and generates the frequency distribution of the column
    described by paramter, most likely EW_r or logN
    '''

    # Open the file
    try:
        d = pd.read_hdf(filename, 'data')
        data = d[paramter]
    except:
        print 'ERROR in bindata-logreq.py while accessing:'
        print '\tFile = {0:s}\n\tParameter = {1:s}'.format(filename, parameter)
        return 1

    numdata = len(data) 
    numbins = int((upperlimit - lowerlimit) / binsize)


    # Bin the data
    h, edges = np.histogram(data, bins=numbins, range=(lowerlimit,upperlimit)) 

    # Determine the errors
    binup, bindn = np.zeros(numbins), np.zeros(numbins)
    for i in range(numbins):

        # For large values, use the Gaussian approximation
        if h[i]>200:
            binup[i] = h[i] + sqrt(h[i])
            bindn[i] = h[i] - sqrt(h[i])
        
        # To be done later
        #elif h[i]>0 and h[i]<=200:
            # Use Possion CDF to computer uncertainties
            #n = h[i]
            #xllim = 0.0
            #xulim = 1000.0*n
            
        else:
            # Placeholder 
            binup[i] = h[i] + sqrt(h[i])
            bindn[i] = h[i] - sqrt(h[i])

    return h, edges, binup, bindn










