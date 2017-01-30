# 
# Description:
#  Determines which cells are the significant contributers 
#  to the EW measurement for all LOS for a given ion/snapshot
#

import numpy as np
import sys
from  math import sqrt, pow, log10
import subprocess as sp
import ew 
import os.path as op
import locate_funcs as lf
import os
import glob
import significant as sg
import pandas as pd
import locate_files as fi


def locateSigCells(galID, expn, ion, ewcut, codeLoc, inc, testing=0):

    """
    Locates the significant cells in a line of sight.
    Remove cells in order of increasing column density as determined
    by los7 until the ew changes by ewcut. These spectra are generated
    with specsynth with no noise. 
    
    Results are reported in a file called
        <galID>.<expn>.<ion>.i<inc>abs_cells.dat
    """

    singleCount = 0

    # Read in the galaxy's box
    boxfile = '../{0:s}_GZa{1:s}.{2:s}.h5'.format(galID,expn,ion)
    box = pd.read_hdf(boxfile, 'data')
    if testing==1:
        print 'Box read in'

    # Read in the LOS info from lines.info
    los_info = np.loadtxt('lines.info',skiprows=2)

    # Determine which transition to use
    # Only need on, use the transition with the weaker oscillator strength since
    # this will retain mroe cells than the other transitions
    wave = fi.transition_name(ion)

    # Open the output file
    outfile = '{0:s}.{1:s}.{2:s}.i{3:d}.abs_cells.h5'.format(galID,expn,ion,int(inc))
    header = ['LOS','D','cellID','redshift','logN','dobbler_b','x', 'y', 'z', 'vx', 'vy', 'vz',
                'galactocentric_d', 'nH', 'temperature', 'cell_size', 'SNII', 'SNIa', 
                'alpha_Zmet', 'ion_density']

    # Write a header to the output file

    numcols = len(header)

    # Create a black row of zeros to build the array with
    d = np.zeros(numcols)

    # Make a version of Mockspec.runpars that has zero SNR
    # Needed for sigcells
    lf.quiet_mockspec()

    # Get a list of LOS that have a sysabs file associated with it
    sysabs_losnum = glob.glob('*los*sysabs')
    sysabs_losnum.sort()
    
    if testing==1:
        print 'Sysabs files aggregated'
        print 'Number of sysabs files: ', len(sysabs_losnum)

    flog = open('sig_cells.log', 'w')
    for i in range(0,len(sysabs_losnum)):

        losnum = sysabs_losnum[i].split('.')[2].split('los')[1]
        num = int(losnum)

        linesfile = '{0:s}.{1:s}.los{2:s}.lines'.format(galID,ion,losnum)

        # Make sure the .lines file has cells in it
        numCells = 0
        with open(linesfile) as f:
            for line in f:
                numCells += 1

        # There is always the galaxy redshift
        if numCells==1:
            # If there are no cells, continue to the next LOS
            continue        

        # Get the impact paramter of this LOS
        imp = los_info[num-1,1]
        
        # Copy the original lines file
        command = 'cp '+linesfile+' '+linesfile+'.tmp'
        sp.call(command, shell=True)

        # Perform the velocity cut
        lf.velcut(linesfile, testing=testing)

        # Find the significant cells
        endCut = sg.sigcells(linesfile,ewcut,codeLoc,flog,wave,testing=testing)

        # Get the properties of the cells
        # Open the lines.final file
        finalLinesFile = linesfile.replace('.lines',
                            '{0:s}.lines.final'.format(wave))
        final_file = open(finalLinesFile)
        final_file.readline()
        for line in final_file:
            l = line.split()
            redshift = float(l[0])
            column = float(l[1])
            doppler = float(l[2])
            cellID = int(float(l[3]))
        
            # Get the cell's properties from the boxfile
            index = cellID-1

            cell_size = box['cell_size'].iloc[index]
            x = box['x'].iloc[index]
            y = box['y'].iloc[index]
            z = box['z'].iloc[index]
            vx = box['vx'].iloc[index]
            vy = box['vy'].iloc[index]
            vz = box['vz'].iloc[index]
            density = box['nH'].iloc[index]
            temperature = box['temperature'].iloc[index]
            snII = box['SNII'].iloc[index]
            snIa = box['SNIa'].iloc[index]
            alphaZ = box['alpha_Zmet'].iloc[index]
            ionDense = box['nIon'].iloc[index]
        
            # Calculate the galactocentric distance
            r = np.sqrt(x*x + y*y + z*z)
            
            # Fill an array of these values
            sig = np.zeros(numcols)
            sig[0] = float(num)
            sig[1] = imp
            sig[2] = cellID
            sig[3] = redshift
            sig[4] = column
            sig[5] = doppler
            sig[6] = x
            sig[7] = y
            sig[8] = z
            sig[9] = vx
            sig[10] = vy
            sig[11] = vz
            sig[12] = r
            sig[13] = density
            sig[14] = temperature
            sig[15] = cell_size
            sig[16] = snII
            sig[17] = snIa
            sig[18] = alphaZ
            sig[19] = ionDense
        
            # Append to the main array
            d = np.vstack((d,sig))

            # Write all to the output file
            #s = '{0:d}'.format(num).ljust(7) +  '{0:.3f}'.format(imp).rjust(7) +  '{0:d}'.format(cellID).rjust(16) + '{0:-.7f}'.format(redshift).rjust(14) + '{0:.3f}'.format(column).rjust(10) + '{0:.3f}'.format(doppler).rjust(13) + '{0:.5e}'.format(r).rjust(20) + '{0:.4f}'.format(density).rjust(12) + '{0:.4f}'.format(temperature).rjust(10) + '{0:.4f}'.format(cell_size).rjust(14) + '{0:.4e}'.format(snII).rjust(19) + '{0:.4e}'.format(snIa).rjust(20) + '{0:.4e}'.format(alphaZ).rjust(17) + '{0:.4e}'.format(ionDense).rjust(17) + '\n'
            #f_out.write(s)

        # Rename the original .lines file back to its full name
        command = 'cp '+linesfile+'.tmp '+linesfile
        sp.call(command, shell=True)
      
    # Write the outfile
    df = pd.DataFrame(d, columns=header)
    df.to_hdf(outfile, 'data', mode='w')

    
    #flog.close()  
    #f_out.close()

    print 'For {0:s}, {1:d} LOS are dominated by one cell'.format(ion, singleCount)





















