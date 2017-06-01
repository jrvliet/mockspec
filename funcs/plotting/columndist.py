
# Plots the column density distribution for significantly absorbing cells

# Usage:
#   python 


from __future__ import print_function
import numpy as np
import pandas as pd
import scipy.optimize as so
import matplotlib.pyplot as plt
import subprocess as sp
import os
import sys

def schecter(l,phi,lstar,alpha):

    schec = (phi/lstar) * pow((l/lstar),alpha) * np.exp(-1.*l/lstar)
    return schec

def column_distribution(ions):

    # Read in the galaxy properties from galaxy.props
    f = open('galaxy.props')
    galID = f.readline().split()[1]
    expn = f.readline().split()[1]
    redshift = f.readline().split()[1]
    mvir = f.readline().split()[1]
    rvir = float(f.readline().split()[1])
    inc = int(float(f.readline().split()[1]))
    f.close()

    fitHeader = 'ion phi lstar alpha sigmaPhi sigmaLstar sigmaAlpha'.split()
    fit = np.zeros((len(ions),len(fitHeader)))
    fit = pd.DataFrame(fit,columns=fitHeader)

    # Command line arguements for binfreq
    column = 19             # Column the data are in (count from 1)
    linear = 0              # Input data type (0=linear, 1=log10)
    binsize = 0.1           # Binsize (equal log10)
    lowerlimit = -3         # Lower bin limit (log10)
    upperlimit = 2          # Upper bin limit (log10)
    header = 1              # Number of header rows

    # Get the location of the code
    pathname = os.path.dirname(sys.argv[0])
    codeLoc = os.path.abspath(pathname)

    fig,axes = plt.subplots(2,2,figsize=(10,10))

    # Loop over ions
    for i,(ion,ax) in enumerate(zip(ions,axes.flatten())):
            
        allfile = './{0:s}/{1:s}.{0:s}.a{2:s}.i{3:d}.ALL.sysabs.h5'.format(ion.name,
                    galID,expn,inc)
        d = pd.read_hdf(allfile,'data')
        tmpfile = 'logNCols.tmp'
        np.savetxt(tmpfile,d['logN'])

        # Run Chris's binning program
        blankCommand = '{0:s}/funcs/plotting/bindata-logfreq {1:s} {2:d} {3:d} {4:f} {5:f} {6:f} {7:d}'
        command = blankCommand.format(codeLoc, allfile, column, linear, binsize, lowerlimit, upperlimit, header)
        command = blankCommand.format(codeLoc, tmpfile, 0, linear, binsize, lowerlimit, upperlimit, header)


        try:
            sp.check_call(command, shell=True)
        except:
            print('Error running binning program in columndist.py')
            break

        # Output will be named <galID>.logfreqbin
        # Rename to <galID>.<expn>.<ion>.ew.logfreqbin
        oldname = 'logNCols.logfreqbin'
        newname = '{0:s}.{1:s}.{2:s}.logN.logfreqbin'.format(galID, expn, ion.name)
        command = 'mv {0:s} {1:s}'.format(oldname, newname)
        try:
            sp.call(command, shell=True)
        except:
            print('Error renaming files in columndist.py')
            break

        # Read in the binned data
        try:
            cols = 'binCenter logF halfBin dlogFdn dlogFup'.split()
            data = pd.read_csv(newname,skiprows=4,sep='\s+',names=cols)
            data = data[data['logF']>-10]
        except IOError:
            print('Error reading in {0:s} in columndist.py'.format(newname))
            raise

        # Fit Schecter to the data
        p0 = [1.,-0.1,-1.0]
        try:
            fitp,fitsig = so.curve_fit(schecter,data['binCenter'],data['logF'],
                                        p0=p0,sigma=data['dlogFup'],absolute_sigma=False)
        except (RuntimeError,TypeError) as e:
            fitp = p0
            figsig = np.identity(3)

        phi,lstar,alpha = fitp
        sigma = np.sqrt(np.diag(fitsig))

        fit['ion'].iloc[i] = ion.name
        fit['phi'].iloc[i] = phi
        fit['lstar'].iloc[i] = lstar
        fit['alpha'].iloc[i] = alpha
        fit['sigmaPhi'].iloc[i] = sigma[0]
        fit['sigmaLstar'].iloc[i] = sigma[1]
        fit['sigmaAlpha'].iloc[i] = sigma[2]
        
        data['fit'] = data['binCenter'].apply(schecter,args=(phi,lstar,alpha))
            
        # Plot the data
        ax.errorbar(data['binCenter'],data['logF'],xerr=data['halfBin'],
                    yerr=[data['dlogFdn'],data['dlogFup']],
                    marker='s',color='k',linestyle='none',label='Data')
        
        # Plot the fit
        ax.plot(data['binCenter'],data['fit'],color='r',label='Fit')

        # Format plot
        ax.set_xlabel('log( Column Density [cm$^{-2}$] )')
        ax.set_ylabel('log ( n() )')
        ax.set_xlim([-3, 2])
        ax.set_ylim([-5, 2])
        ax.legend(frameon=False,loc='lower left')

    fig.tight_layout()
    fig.subplots_adjust(top=0.92)
    fig.suptitle(r'{0:s}, a={1:s}, Rvir={2:.1f} kpc'.format(galID, expn, rvir))
    s = '{0:s}_{1:s}_columndist.pdf'.format(galID, expn)
    fig.savefig(s, bbox_inches='tight')
    plt.close(fig)

    fitName = '{0:s}.a{1:s}.i{2:d}.logNFits.dat'.format(galID,expn,inc)
    fit.to_csv(fitName,index=False)





