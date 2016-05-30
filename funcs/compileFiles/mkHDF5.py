
'''
Codes to convert the ALL.sysabs files and 
abs_cells files to HDF5 format for
compactness
'''

import numpy as np
import pandas as pd
import glob

def sysbas_to_hdf5(codeLoc):

    '''
    Converts the ALL.sysabs file into the HDF5 format
    '''
    
    # Get the name of the file 
    allfile = glob.glob('*.ALL.sysabs')
    af = allfile.split('.')
    galID = af[0]
    ion = af[1]
    expn = af[2].split('a')[1]
    inc = af[3].split('i')[1]

    # Create the name of the HDF5 file
    hdf5file = allfile + '.hdf5'

    # Open the file in append mode
    hdf = HDFStore(hdf5file)



        






