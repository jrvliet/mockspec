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


def locateSigCells(galID, expn, ion, ewcut, codeLoc, inc, testing=0):

    """
    Locates the significant cells in a line of sight.
    Remove cells in order of increasing column density as determined
    by los7 until the ew changes by ewcut. These spectra are generated
    with specsynth with no noise. 
    
    Results are reported in a file called
        <galID>.<expn>.<ion>.abs_cells.dat
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
    out_file = '{0:s}.{1:s}.{2:s}.abs_cells.dat'.format(galID,expn,ion)
    #out_file = galID+'.'+expn+'.'+ion+'.abs_cells.dat'
    f_out = open(out_file, 'w')

    # Write a header to the output file
    header = 'LOS \t Imp Param      Cell ID   Redshift        logN'\
            '    Doppler b    Galactocentric d      log nH     log T'\
            '     Cell Size     SNII mass frac      SNIa mass frac'\
            '       alpha_Zmet     Ion Density\n'
    f_out.write(header)


    # Make a version of Mockspec.runpars that has zero SNR
    # Needed for sigcells
    lf.quiet_mockspec()

    # Get a list of LOS that have a sysabs file associated with it
    sysabs_losnum = glob.glob('*los*sysabs')
    sysabs_losnum.sort()
    
#    sysabs_losnum = []
#    for i in range(0,1000):
#        losnum = str(i).zfill(4)
#        filename = galID+'.'+ion+'.los'+losnum+'.sysabs'
        
        # Check to see if the file exists
#        if op.isfile(filename):
#            sysabs_losnum.append(losnum)

    if testing==1:
        print 'Sysabs files aggregated'
        print 'Number of sysabs files: ', len(sysabs_losnum)

    for i in range(0,len(sysabs_losnum)):

        losnum = sysabs_losnum[i].split('.')[2].split('los')[1]

        if testing==1:
            print 'LOS num: ', losnum
            if int(losnum)!=234:
                continue
        
        num = int(losnum)
        linesfile = galID+'.'+ion+'.los'+losnum+'.lines'

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
        if testing==1:
            print '\tPerforming velocity cut'
        lf.velcut(linesfile, testing=testing)

        # Find the significant cells
        print linesfile
        if testing==1:
            print '\t Finding significant cells'
        endCut = sg.sigcells(linesfile, ewcut, codeLoc, testing=testing)
#        singleCount += lf.sigcells(linesfile, ewcut, codeLoc, testing=testing)
        print linesfile
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
            density = np.log10(box[index,7])
            temperature = np.log10(box[index,8])
            snII = box[index,9]
            snIa = box[index,10]
            alphaZ = box[index,16]
            ionDense = box[index,13]

            # Calculate the galactocentric distance
            r = np.sqrt(x*x + y*y + z*z)
        
            # Write all to the output file
            s = '{0:d}'.format(num).ljust(7) +  '{0:.3f}'.format(imp).rjust(7) +  '{0:d}'.format(cellID).rjust(16) + '{0:-.7f}'.format(redshift).rjust(14) + '{0:.3f}'.format(column).rjust(10) + '{0:.3f}'.format(doppler).rjust(13) + '{0:.5e}'.format(r).rjust(20) + '{0:.4f}'.format(density).rjust(12) + '{0:.4f}'.format(temperature).rjust(10) + '{0:.4f}'.format(cell_size).rjust(14) + '{0:.4e}'.format(snII).rjust(19) + '{0:.4e}'.format(snIa).rjust(20) + '{0:.4e}'.format(alphaZ).rjust(17) + '{0:.4e}'.format(ionDense).rjust(17) + '\n'
            f_out.write(s)

        # Rename the original .linse file back to its full name
        command = 'cp '+linesfile+'.tmp '+linesfile
        sp.call(command, shell=True)
        
    f_out.close()

    print 'For {0:s}, {1:d} LOS are dominated by one cell'.format(ion, singleCount)

