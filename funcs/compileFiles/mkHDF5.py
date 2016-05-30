
'''
Codes to convert the ALL.sysabs files and 
abs_cells files to HDF5 format for
compactness
'''

import numpy as np
import pandas as pd
import glob

def sysabs_to_hdf5(codeLoc):

    '''
    Converts the ALL.sysabs file into the HDF5 format
    '''
    
    # Get the name of the file 
    allfile = glob.glob('*.ALL.sysabs')

    # Create the name of the HDF5 file
    hdf5file = allfile + '.hdf5'

    # Get the header
    with open(allfile) as f:
        header = f.readline().strip().split()

    # Read in the data
    data = np.loadtxt(allfile, skiprows=1)
    
    # WRite data to HDF file
    df = pd.DataFrame(data, columns=header)
    df.to_hdf(hdffile, 'data', mode='w')


def abscells_to_hdf5(codeLoc):

    '''
    Converts the abs_cells.dat file into the HDF5 format
    '''
    
    # Get the name of the file 
    filename = glob.glob('*.abs_cells.dat')

    # Create the name of the HDF5 file
    hdf5file = filename.replace('dat','hdf5')

    # Get the header
    with open(filename) as f:
        header = f.readline().strip().split()

    # Read in the data
    data = np.loadtxt(filename, skiprows=1)
    
    # WRite data to HDF file
    df = pd.DataFrame(data, columns=header)
    df.to_hdf(hdffile, 'data', mode='w')


def gasbox_to_hdf5(codeLoc):

    '''
    Converts the ion boxes file into the HDF5 format
    '''
    
    # Get the name of the file 
    files = glob.glob('*GZ*.txt')

    for filename in files:

        # Create the name of the HDF5 file
        hdf5file = filename.replace('txt','hdf5')

        # Get the header
        with open(filename) as f:
            header = f.readline().strip().split()

        # Read in the data
        data = np.loadtxt(filename, skiprows=1)
        
        # WRite data to HDF file
        df = pd.DataFrame(data, columns=header)
        df.to_hdf(hdffile, 'data', mode='w')






        






