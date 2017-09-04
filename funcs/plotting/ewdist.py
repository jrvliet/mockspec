
# Plots the EW distribution 

# Usage:
#   python 


from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('agg')
import pandas as pd
import subprocess as sp
import os
import sys
import scipy.optimize as so


def schechter(l, phi, lstar, alpha):

    schec = (phi/lstar) * pow((l/lstar), alpha) * np.exp(-1.0*l/lstar)
    return schec    

#def ew_distribution(ions):
def ewdist(ions):

    # Initial guess for Schechter parameters
    # Phi = 1
    # Lstar = -0.100
    # Alpha = -1.0
    p0 = [1.000, -0.100, -1.0]

    # Read in the galaxy properties from galaxy.props
    f = open('galaxy.props')
    galID = f.readline().split()[1]
    expn = f.readline().split()[1]
    redshift = f.readline().split()[1]
    mvir = f.readline().split()[1]
    rvir = float(f.readline().split()[1])
    inc = int(float(f.readline().split()[1]))
    f.close()

    # Open output file
    outfile = 'ewdist_schechter_{0:s}_a{1:s}_i{2:d}.out'.format(galID, expn, inc)
    fout = open(outfile, 'w')
    header = 'Ion\tPhi\t\t\tW*\t\t\tAlpha\n'
    fout.write(header)

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

    fig,axes = plt.subplots(2,2,figsize=(10,10))
    # Loop over ions
    for ion,ax in zip(ions,axes.flatten()):
        
        allfile = './{0:s}/{1:s}.{0:s}.a{2:s}.i{3:d}.ALL.sysabs.h5'.format(
                    ion.name,galID,expn,inc)
        
        d = pd.read_hdf(allfile, 'data')
        ew = d['EW_r']
        tmpfile = 'ewCols.tmp'
        np.savetxt(tmpfile, ew)

        # Run Chris's binning program
        blankCommand = '{0:s}/funcs/plotting/bindata-logfreq {1:s} {2:d} {3:d} {4:f} {5:f} {6:f} {7:d}'
        command = blankCommand.format(codeLoc, tmpfile, 0, linear, binsize, 
                                        lowerlimit, upperlimit, header)
        sp.call(command, shell=True)

        # Output will be named <galID>.logfreqbin
        # Rename to <galID>.<expn>.<ion>.ew.logfreqbin
        oldname = '.logfreqbin'.format(galID)
        oldname = 'ewCols.logfreqbin'
        newname = '{0:s}.{1:s}.{2:s}.ew.logfreqbin'.format(galID, expn, ion.name)
        command = 'mv {0:s} {1:s}'.format(oldname, newname)
        sp.call(command, shell=True)

        # Read in the binned data
        try:
            #data = np.loadtxt(newname, skiprows=4)
            cols = 'binCenter logF halfBin dlogFdn dlogFup'.split()
            data = pd.read_csv(newname,skiprows=4,sep='\s+',names=cols)
            data = data[data['logF']>-10]
        except IOError:
            print('Error in ew_distribution in ewdist while reading {0:s}'.format(newname))
            raise

        # Fit a Schecter function to the data
        try:
            fitp,fitsig = so.curve_fit(schechter,data['binCenter'],data['logF'],p0=p0,
                                    sigma=data['dlogFup'],absolute_sigma=False)
        except RuntimeError:
            fitp = p0
            fitsig = np.identity(3)
    
        phi,lstar,alpha = fitp
        sigma = np.sqrt(np.diag(fitsig))

        data['fit'] = data['binCenter'].apply(schechter,args=(fitp[0],fitp[1],fitp[2]))

        # Write results to file
        s = '{0:s}\t{1:f} +/- {2:f}\t{3:f} +/- {4:f}\t{5:f} +/- {6:f}\n'.format(
            ion.name, phi, sigma[0], lstar, sigma[1], alpha, sigma[2])
        fout.write(s)

        # Plot the data
        ax.errorbar(data['binCenter'],data['logF'],xerr=data['halfBin'],
                    yerr=[data['dlogFdn'],data['dlogFup']], 
                    marker='s',color='k',linestyle='none',label='Data')
        
        # Plot the fit
        ax.plot(data['binCenter'],data['fit'],color='r',label='Fit')

        # Format plot
        ax.set_xlabel('log( EW [$\AA$] )')
        ax.set_ylabel('log ( n(EW) )')
        ax.set_xlim([-3, 2])
        ax.set_ylim([-5, 2])
        ax.legend(frameon=False, loc='lower left', prop={'size':8})

        sp.call('rm {0:s}'.format(tmpfile), shell=True)

#        xdata, ydata, err = [], [], []
#        for i in range(0,len(freq)):
#            if freq[i]>-50:
#                x = pow(10.0,binCenter[i])
#                y = pow(10.0,freq[i])
#                xdata.append(x)
#                ydata.append(y)
#                errorU = pow(10.0,errUp[i])
#                errorD = pow(10.0,errDown[i])
#                meanError = (errorU+errorD)/2.0
#                err.append( meanError )
#
#        (phi, lstar, alpha), paramCovariance = curve_fit(schechter, xdata, ydata, p0=p0, 
#                                                         sigma=err, absolute_sigma=True)
#
#
#        # Plot the data
#        xerrbin = pow(10.0,halfbin)
#        yerrbinDown = pow(10.0, errDown)
#        yerrbinUp = pow(10.0, errUp)
#        ax.errorbar(binCenter, freq, xerr=halfbin, yerr=[errDown,errUp], 
#                    marker='s',color='k',
#                    linestyle='none', label='Data')
#
#        # Overplot fit
#        y = []
#        x = []
#        for l in xdata:
#            val = schechter(l, phi, lstar, alpha)
#            y.append(np.log10(val))
#            x.append(np.log10(l))
#            #x.append(l)
#            #y.append(val)

    fig.tight_layout()
    fig.subplots_adjust(top=0.92)
    fig.suptitle(r'{0:s}, a={1:s}, Rvir={2:.1f} kpc, i={3:d}'.format(galID, expn, rvir, inc))
    s = '{0:s}_a{1:s}_i{2:d}_ewdist.pdf'.format(galID, expn, inc)
    fig.savefig(s, bbox_inches='tight')
    fout.close()





