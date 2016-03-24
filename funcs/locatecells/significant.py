
import numpy as np
import subprocess as sp
import locate_funcs as lf
import locate_files as fi
from ew import findEW
import sys


def search(start, end, ewcut, lines_z, lines_b, lines_N, lines_ID, redshift, 
            specCommand, linesfile, specfile, negVelLimit, posVelLimit, 
            ewVelcut, flog, depth):

    '''
    Recursively search for significant cells
    '''
    maxDepth = 500   

    # Check if are at the end
    if depth>maxDepth:
        flog.write('{0:d}\t{1:d}\tMax Depth Reached\n'.format(start,end))
        return end
    if abs(start-end)<=1:
        flog.write('{0:d}\t{1:d}\tDone\n'.format(start,end))
        return end
    else:

        depth += 1

        # Get midpoint
        mid = int((end-start)/2)
        index = mid+start
        mid = index
        
        # Cut the lines
        cut_z = lines_z[:index]
        cut_N = lines_N[:index]
        cut_b = lines_b[:index]
        cut_ID = lines_ID[:index]

        # Write the new array to the lines file
        s = '{0:>8.7f}\t{1:>8f}\t{2:>8f}\t{3:>8d}\n'
        with open(linesfile, 'w') as f:            
            f.write('{0:.16f}\n'.format(redshift))
            for i in range(0,len(cut_z)):
                f.write(s.format(cut_z[i], cut_N[i],
                                 cut_b[i], int(cut_ID[i])))
        # Run specsynth
        sp.call(specCommand, shell=True)

        # Find the new EW of the cut down list of cells
        wavelength, velocity, flux = np.loadtxt(specfile, usecols=(0,1,2),
                                                unpack=True)
        ew = findEW(wavelength, velocity, flux, negVelLimit, posVelLimit)
        ewdiff = abs( (ewVelcut - ew) / ewVelcut)*100
        flog.write('{0:d}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\t{5:.3f}\t{6:.3f}\n'.format(len(lines_z), start, end, mid, len(cut_z), ew, ewdiff))
    
        if ewdiff<ewcut:
            # Not deep enough cut
            # Only send in top half
            end = mid
            search(start, end, ewcut, lines_z, lines_b, lines_N, lines_ID, 
                    redshift, specCommand, linesfile, specfile, negVelLimit, 
                    posVelLimit, ewVelcut, flog, depth)
        else:
            # Cut is too deep
            # Need to include more of cut cells 
            start = mid
            search(start, end, ewcut, lines_z, lines_b, lines_N, lines_ID, 
                    redshift, specCommand, linesfile, specfile, negVelLimit, 
                    posVelLimit, ewVelcut, flog, depth)
    return end

def sigcells(linesfile, ewcut, codeLoc, flog, testing=0):
    '''
    Determines which cells are the significiant contributers
    to the EW measurement. Cells are removed until the EW of 
    the absorption feature (with zero SNR) is different from the 
    EW of the absorption (0 SNR) from the full, original list of 
    cells by a percentage deteremined by ewcut
    '''

    singleCellCount = 0       # Counts number of LOS dominated by a single cell
    
    # Get the info from the filename
    galID  = linesfile.split('.')[0]
    ion    = linesfile.split('.')[1]
    losnum = (linesfile.split('.')[2]).split('los')[1]

    flog.write('\n{0:s}\n'.format(losnum))

    # Read in the linesfile
    with open(linesfile+'.velcut','r') as f:
        redshift = float(f.readline().strip())
    cell_z, cell_N, cell_b, cell_ID = np.loadtxt(linesfile+'.velcut', skiprows=1,
                                    usecols=(0,1,2,3), unpack=True)

    # Sort the arrays in decreasing column density
    # as long as there are more than one cell
    if type(cell_z) is np.ndarray:
        inds = np.argsort(cell_N)[::-1]
        cell_z = np.array(cell_z)[inds]
        cell_N = np.array(cell_N)[inds]
        cell_b = np.array(cell_b)[inds]
        cell_ID = np.array(cell_ID)[inds]

    if testing==1:
        print 'In sigcells, number of velcut cells read in: ', len(cell_z)

    # Get the EW from sysabs
    negVelLimit, posVelLimit, ewSysabs = lf.vel_limits(linesfile)

    # Create a los.list file containing only this LOS
    with open('los_single.list', 'w') as fLos:
        datfile = linesfile.replace('lines', 'dat') + '\n'
        fLos.write(datfile)
    
    # Rename the velcut .lines to remove velcut from name, 
    # so it will be used by specsynth
    command = 'cp '+linesfile+'.velcut '+linesfile
    sp.call(command, shell=True)

    # Run specsynth on the velcut lines list
    specCommand = codeLoc+'/funcs/mkspec/specsynth los_single.list Mockspec_0SNR.runpars'
    sp.call(specCommand, shell=True)

    # Get the EW of this noise-less spectra
    bluewave, redwave = fi.transition_name(ion, codeLoc)
    specFileBase = '{0:s}.{1:s}.los{2:s}.{3:s}.spec'
    specfile = specFileBase.format(galID, ion, losnum, bluewave)
    redspecfile = specFileBase.format(galID, ion, losnum, redwave)

    # Copy the initial quiet spectra
    command = 'cp {0:s} {0:s}.velcutclean'.format(specfile)
    sp.call(command, shell=True)
    command = 'cp {0:s} {0:s}.velcutclean'.format(redspecfile)
    sp.call(command, shell=True)
    
    # Read in the spectra data
    specdata = np.loadtxt(specfile)
    wavelength = specdata[:,0]
    velocity = specdata[:,1]
    flux = specdata[:,2]
    
    ew = findEW(wavelength, velocity, flux, negVelLimit, posVelLimit)
    ewdiff = abs( (ewSysabs - ew) / ewSysabs )
    ew_velcut_lines = ew
    ewVelcut = ew

    #################################################################
    #                                                               #
    #     Remove cells until the EW differs from full .lines        #
    #     value by more than ewcut                                  #
    #                                                               #
    #################################################################
    
    # See how many cells are in the velcut file
    with open(linesfile+'.velcut') as fcut:
        for i,l in enumerate(fcut):
            pass
    numcells = i   # One less than the actual number of lines due to the
                   # absorption redshift in the first line
    
    # If there is only one cell in the file, stop. Clearly that cell
    # is responsible for all the absorption
    s = '{0:>8.7f}\t{1:>8f}\t{2:>8f}\t{3:>8d}\n'
    with open(linesfile+'.final', 'w') as f:
        f.write('{0:.16f}\n'.format(redshift))

        if numcells == 1:
            singleCellCount += 1
            sigEnd = 1    
            f.write(s.format(cell_z, cell_N, cell_b, int(cell_ID)))

        else:
            start = 0
            end = len(cell_z)
            depth = 0
            sigEnd = search(start, end, ewcut, cell_z, cell_b, cell_N, cell_ID, 
                redshift, specCommand, linesfile, specfile, negVelLimit, 
                posVelLimit, ewVelcut, flog, depth)

            for i in range(0,sigEnd):
                f.write(s.format(cell_z[i], cell_N[i], cell_b[i], int(cell_ID[i])))
            



