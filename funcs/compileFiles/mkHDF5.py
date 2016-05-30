
'''
Codes to convert the ALL.sysabs files and 
abs_cells files to HDF5 format for
compactness
'''

import numpy as np
import pandas as pd
import glob
import losSummary as ls

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
        data = np.loadtxt(filename, skiprows=2)
        
        # WRite data to HDF file
        df = pd.DataFrame(data, columns=header)
        df.to_hdf(hdf5file, 'data', mode='w')





def genSummaries(galID, expn, incline, ions, numlos):

    '''
    A controlling function to generate LOS summary files.
    Uses the functions in losSummary.py

    Output file has the name:
        <galID>_a<expn>_i<incline>_<ion>_cellSummary.txt
    '''

    # Get the IDs of cells along the LOS
    probbedIDsAll = []
    for i in range(numlos):
        losnum = str(i+1).zfill(4)
        probbedIDsAll.append(ls.num_along_los(losnum))

    # Read in the impact parameters and azimuthal angles
    imp, phi = np.loadtxt('lines.info', skiprows=2,
                        usecols=(1,2), unpack=True)
    
    for ion in ions:

        outfile = '{0:s}_a{1:s}_i{2:d}_{3:s}_cellSummary.txt'.format(
                    galID,expn,incline,ion)

        # Write the header
        header = ('LOS\tImpact\tPhi\tProbbed\tProbbed_in_halos\t'
                    'In_Lines\tLines_in_halos\tSignificant\t'   
                    'Sig_in_halos\n')
        f = open(outfile, 'w')
        f.write(header)
                
        s = ('{0:d}\t{1:.3f}\t{2:.3f}\t{3:d}\t{4:d}\t{5:d}\t{6:d}\t'
                '{7:d}\t{8:d}\n')
    
        # Loop over the lines of sight
        for i in range(numlos):
            
            probbedIDs = probbedIDsAll[i]
            losnum = str(i+1).zfill(4)

            # Get the id of cells in the .lines file for this LOS 
            linesIDs = num_in_lines(galID, ion, losnum)

            # Get the id of the significant cells along this LOS
            sigIDs = num_significant(galID, expn, ion, incline, losnum)

            # Determine the number of these cells that are in subhalos
            probbedinSub, linesinSub, siginSub = ls.num_in_subhalos(galID,
                expn, ion, incline, losnum, probbedCells, linesCells, 
                sigCells)


            # Format the output string
            f.write(s.format(i+1, imp[i], phi[i], len(probbedIDs), 
                probbedinSub, len(linesIDs), linesinSub, len(sigIDs),
                siginSub))

        f.close()








