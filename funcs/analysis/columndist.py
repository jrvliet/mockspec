
# Plots the column density distribution for significantly absorbing cells

# Usage:
#   python 


import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
import os
import sys

def column_distribution(ions)

    # Read in the galaxy properties from galaxy.props
    f = open('galaxy.props')
    galID = f.readline().split()[1]
    expn = f.readline().split()[1]
    redshift = f.readline().split()[1]
    mvir = f.readline().split()[1]
    rvir = float(f.readline().split()[1])
    f.close()


    # Command line arguements for binfreq
    column = 6              # Column the data are in (count from 1)
    linear = 0              # Input data type (0=linear, 1=log10)
    binsize = 0.1           # Binsize (equal log10)
    lowerlimit = -3         # Lower bin limit (log10)
    upperlimit = 2          # Upper bin limit (log10)
    header = 1              # Number of header rows

    # Get the location of the code
    pathname = os.path.dirname(sys.argv[0])
    codeLoc = os.path.abspath(pathname)

    # Loop over ions
    for ion in ions:
            
        allfile = './'+ion+'/'+galID+'.'+ion+'.a'+expn+'.ALL.sysabs'
        
        
        absfile = './{0:s}/{1:s}.{2:s}.a{0:s}.abs_cells.dat'.format(ion, galID, expn)

        # Run Chris's binning program
        blankCommand = '{0:s}/bindata-logfreq {1:s} {2:d} {3:d} {4:f} {5:f} {6:f} {7:d}'
        command = blankCommand.format(codeLoc, absfile, column, linear, binsize, lowerlimit, upperlimit, header)

        sp.call(command, shell=True)

        # Output will be named <galID>.logfreqbin
        # Rename to <galID>.<expn>.<ion>.ew.logfreqbin
        oldname = '.logfreqbin'.format(galID)
        newname = '{0:s}.{1:s}.{2:s}.ew.logfreqbin'.format(galID, expn, ion)
        command = 'mv {0:s} {1:s}'.format(oldname, newname)
        sp.call(command, shell=True)

        # Read in the binned data
        data = np.loadtxt(newname, skiprows=4)
        binCenter = data[:,0]
        freq = data[:,1]
        halfbin = data[:,2]
        errDown = -1.0*data[:,3]
        errUp   = data[:,4]

        # Plot the data
        subplotnum = 221+ion_list.index(ion)
        plt.subplot(subplotnum)
        plt.errorbar(binCenter, freq, xerr=halfbin, yerr=[errDown,errUp], linestyle='none')
        plt.xlabel('log( Column Density [cm$^{-2}$] )')
        plt.ylabel('log ( n() )')
        plt.xlim([-3, 2])
        plt.ylim([-5, 2])

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.suptitle(r'{0:s}, a={1:s}, Rvir={2:.1f} kpc'.format(galID, expn, rvir))
    s = '{0:s}_{1:s}_columndist.pdf'.format(galID, expn)
    plt.savefig(s, bbox_inches='tight')






