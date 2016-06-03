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
    boxfile = galID+'_GZa'+expn+'.'+ion+'.txt'
    box = np.loadtxt(boxfile, skiprows=2)
    if testing==1:
        print 'Box read in'

    # Read in the LOS info from lines.info
    los_info = np.loadtxt('lines.info',skiprows=2)

    # Open the output file
    out_file = '{0:s}.{1:s}.{2:s}.i{3:d}.abs_cells.h5'.format(galID,expn,ion,int(inc))
    #f_out = open(out_file, 'w')
    header = ['LOS','D','cellID','redshift','logN','dobbler_b','x', 'y', 'z', 'vx', 'vy', 'vz',
                'galactocentric_d', 'log_nH', 'log_T', 'cell_size', 'SNII', 'SNIa', 'alpha_Zmet', 'ion_density']

    # Write a header to the output file

    numcols = len(header)

    # Create a black row of zeros to build the array with
    d = np.zeros(numcols)

    #header = 'LOS \t Imp_Param      CellID    Redshift        logN'\
    #        '    Doppler_b    Galactocentric_d      log_nH     log_T'\
    #        '     Cell_Size     SNII_mass_frac      SNIa_mass_frac'\
    #        '       alpha_Zmet     Ion_Density\n'
    #f_out.write(header)


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
        endCut = sg.sigcells(linesfile, ewcut, codeLoc, flog, testing=testing)
        # Get the properties of the cells
        # Open the lines.final file
        final_file = open(linesfile+'.final')
        final_file.readline()
        for line in final_file:
            l = line.split()
            redshift = float(l[0])
            column = float(l[1])
            doppler = float(l[2])
            cellID = int(float(l[3]))
        
            # Get the cell's properties from the boxfile
            index = cellID-1
            cell_size = box[index,0]
            x = box[index,1]
            y = box[index,2]
            z = box[index,3]
            vx = box[index,4]
            vy = box[index,5]
            vz = box[index,6]
            density = np.log10(box[index,7])
            temperature = np.log10(box[index,8])
            snII = box[index,9]
            snIa = box[index,10]
            alphaZ = box[index,16]
            ionDense = box[index,13]

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
            d = np.vstack(d,sig)

            # Write all to the output file
            #s = '{0:d}'.format(num).ljust(7) +  '{0:.3f}'.format(imp).rjust(7) +  '{0:d}'.format(cellID).rjust(16) + '{0:-.7f}'.format(redshift).rjust(14) + '{0:.3f}'.format(column).rjust(10) + '{0:.3f}'.format(doppler).rjust(13) + '{0:.5e}'.format(r).rjust(20) + '{0:.4f}'.format(density).rjust(12) + '{0:.4f}'.format(temperature).rjust(10) + '{0:.4f}'.format(cell_size).rjust(14) + '{0:.4e}'.format(snII).rjust(19) + '{0:.4e}'.format(snIa).rjust(20) + '{0:.4e}'.format(alphaZ).rjust(17) + '{0:.4e}'.format(ionDense).rjust(17) + '\n'
            #f_out.write(s)

        # Rename the original .lines file back to its full name
        command = 'cp '+linesfile+'.tmp '+linesfile
        sp.call(command, shell=True)
      
    # Write the outfile
    df = pdf.DataFrame(d, columns=header)
    df.to_hdr(outfile, 'data', mode='w')

    
    #flog.close()  
    #f_out.close()

    print 'For {0:s}, {1:d} LOS are dominated by one cell'.format(ion, singleCount)





















