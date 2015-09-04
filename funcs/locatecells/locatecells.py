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


def locateSigCells(galID, expn, ion, ewcut, codeLoc):

    """
    Locates the significant cells in a line of sight.
    Remove cells in order of increasing column density as determined
    by los7 until the ew changes by ewcut. These spectra are generated
    with specsynth with no noise. 
    
    Results are reported in a file called
        <galID>.<expn>.<ion>.abs_cells.dat
    """

    
    # Read in the galaxy's box
    boxfile = galID+'_GZa'+expn+'.'+ion+'.txt'
    box = np.loadtxt(boxfile, skiprows=2)

    # Read in the LOS info from lines.info
    los_info = np.loadtxt('lines.info',skiprows=2)

    # Open the output file
    out_file = galID+'.'+expn+'.'+ion+'.abs_cells.dat'
    f_out = open(out_file, 'w')

    # Write a header to the output file
    header = 'LOS \t Imp Param      Cell ID   Redshift        logN'\
            '    Doppler b    Galactocentric d      log nH     log T'\
            '     Cell Size     SNII mass frac      SNIa mass frac'\
            '       alpha_Zmet\n'
    f_out.write(header)



    # Get a list of LOS that have a sysabs file associated with it
    sysabs_losnum = []
    for i in range(0,1000):
        
        losnum = str(i).zfill(4)
        filename = galID+'.'+ion+'.los'+losnum+'.sysabs'

        # Check to see if the file exists
        if op.isfile(filename):
            sysabs_losnum.append(losnum)

    for i in range(0,len(sysabs_losnum)):

        losnum = sysabs_losnum[i]
        num = int(losnum)
        linesfile = galID+'.'+ion+'.los'+losnum+'.lines'
        
        # Get the impact paramter of this LOS
        imp = los_info[num-1,1]
        
        # Copy the original lines file
        command = 'cp '+linesfile+' '+linesfile+'.tmp'
        sp.call(command, shell=True)

        # Perform the velocity cut
        lf.velcut(linesfile)

        # Find the significant cells
        lf.sigcells(linesfile, ewcut, codeLoc)

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

            # Calculate the galactocentric distance
            r = np.sqrt(x*x + y*y + z*z)
        
            # Write all to the output file
            s = '{0:d}'.format(num).ljust(7) +  '{0:.3f}'.format(imp).rjust(7) +  '{0:d}'.format(cellID).rjust(16) + '{0:-.7f}'.format(redshift).rjust(14) + '{0:.3f}'.format(column).rjust(10) + '{0:.3f}'.format(doppler).rjust(13) + '{0:.5e}'.format(r).rjust(20) + '{0:.4f}'.format(density).rjust(12) + '{0:.4f}'.format(temperature).rjust(10) + '{0:.4f}'.format(cell_size).rjust(14) + '{0:.4e}'.format(snII).rjust(19) + '{0:.4e}'.format(snIa).rjust(20) + '{0:.4e}'.format(alphaZ).rjust(17) + '\n'

            f_out.write(s)

        # Rename the original .linse file back to its full name
        command = 'cp '+linesfile+'.tmp '+linesfile
        sp.call(command, shell=True)
        
    f_out.close()



