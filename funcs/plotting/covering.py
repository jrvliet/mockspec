#!/usr/bin/python
#
# Filename: coverfrac_master.py
# Version: 1
#
# Author: Jacob Vander Vliet
# Date: 27/07/2014
#
# Compute the non-cumulative covering fraction
# Cut the impact space into annulli
# Add a loop over all EW cuts
# Add in errors
# Fixed the name of galaxies
# Works on the master ALL files that combine results across 
# multiple snapshots
#
# Usage:
#  python generic_covering <galID> <expn> <rvir>
#
# Bulk flag: if=1, will plot all four ions on one plot


from __future__ import print_function
import numpy as np
import pandas as pd
import scipy.special as sc
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('agg')
import sys
import math
import matplotlib
import subprocess as sp
import os

# Function for determining vertical error bars
# Uses incomplete Beta functions
def funcup(x,n1,n2,cl):
    a = n1 + 1.0
    b = n2
    return cl - sc.betainc(a,b,x)

def funcdn(x,n1,n2,cl):
    a = n2 + 1.0
    b = n1
    return cl - sc.betainc(a,b,x)


def covering_fraction(ions):

    legsize = 10

    line=['-','.']

    # Read in galaxy properties from galaxy.props
    f = open('galaxy.props')
    galID = f.readline().split()[1]
    expn = f.readline().split()[1]
    redshift = f.readline().split()[1]
    mass = f.readline().split()[1]
    rad = f.readline().split()[1]
    incline = f.readline().split()[1]

    mvir = float(mass)
    rvir = float(rad)
    inc = int(float(incline))

    CL = 0.8413
    tol = 1e-4
    pmin = 0.0
    pmax = 1.0
    iter = 1000

    # Define the EW cuts to use
    #ewcut = float(sys.argv[1])
    ewcut = 0.0

    # Define the radial bins to use
    nbins = 15
    binsize = 1.5/nbins

    # Determine the bin edges
    Dmax = []
    Dmin = [0.0]
    for j in range(0,nbins):
        Dmax.append( (j+1)*binsize )
        Dmin.append( (j+1)*binsize )
    del Dmin[-1]


    for ion in ions:
        i = 0
        allfile = './{0:s}/{1:s}.{0:s}.a{2:s}.i{3:d}.ALL.sysabs.h5'.format(ion.name, 
                                                                   galID, expn, inc)
        absimpact=[]
        covering=[]
        imp=[]
        horizerrneg=[]
        horizerrpos=[]
        verterrneg=[]
        verterrpos=[]

        # Read in the allfile
        print(os.getcwd())
        try:
            #data = np.loadtxt(allfile, skiprows=1)
            data = pd.read_hdf(allfile,'data')
        except IOError: 
            print('Error in covering.py reading in {0:s}'.format(allfile))
            raise 

        # Loop over the different impact parameters
        for j in range(0,len(Dmax)):
            maxrad = Dmax[j]
            minrad = Dmin[j]

            hit = 0.0
            total = 0.0
            xaxispoint = 0.0
            zerobinabs = []
            zerobin = []
            
            # Loop over all lines from the .sysabs file to get the
            # number of lines with significant absorption (SL > 3)
            # with an impact parameter between minrad and maxrad
            for k in range(0,len(data)):
                impact = data['D'].iloc[k] / rvir
                absimpact.append(impact)
                width = data['EW_r'].iloc[k]
                
                if impact>minrad and impact<maxrad:
                    if width>ewcut:
                        hit += 1
                        xaxispoint += impact
                        total+=1
                        if minrad==0.0:
                            zerobinabs.append( data['los'].iloc[k] )
                    else:
                        total+=1
            
            fraction = hit/total
            covering.append(fraction)
            # Determine the x location of the point. 
            # If hit>0, then the location is the
            # average impact parameter in that bin. 
            # If hit=0, then the location is the
            # midpoint of the bin.
            if hit > 0.0:
                imp.append(xaxispoint/hit)  # location of point on x-axis
            else:
                imp.append((maxrad-minrad)/2.0 + minrad) 
            
            # Determine error bars
            #   Horizontal error bars are width of bin
            horizerrneg.append(imp[j]-minrad)
            horizerrpos.append(maxrad-imp[j])
        
            #   Vertical error bars are found using the incomplete beta function
            top = opt.brentq(funcup,0.0,1.0,args=(hit,total-hit,CL),xtol=tol,maxiter=iter)
            bot = 1.0 - opt.brentq(funcdn,0.0,1.0,args=(hit,total-hit,CL),xtol=tol,maxiter=iter)

            if fraction == 1.0:
                top = 1.0
            elif fraction == 0.0:
                bot = 0.0

            verterrpos.append(top-fraction)
            verterrneg.append(fraction-bot)

            
        subplotnum = 221+ions.index(ion)
        plt.subplot(subplotnum)
        plt.errorbar(imp,covering,xerr=[horizerrneg,horizerrpos],
                     yerr=[verterrneg,verterrpos], linestyle='none')
        i+=1
        plt.xlabel('D / R$_{vir}$')
        ylab = 'C$_f$ ['+ion.name+']'
        plt.ylabel(ylab)
        plt.ylim([-0.1,1.1])
        plt.xlim([min(Dmin),1.5])
        #plt.gca().get_frame().set_linewidth(2) 
        
    #plt.subplots_adjust(wspace=0.3, hspace=0.3)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.suptitle(r'{0:s}, a={1:s}, Rvir={2:.1f} kpc, i={3:d}, EWcut={4:.1f}$\AA$'.format(galID, expn, rvir, inc, ewcut))
    name = '{0:s}_a{1:s}_i{2:d}_{3:.1f}mA_covering.pdf'.format(galID, expn, inc, ewcut*1000)
    plt.savefig(name, bbox_inches='tight')





