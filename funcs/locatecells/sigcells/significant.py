
from __future__ import print_function
import numpy as np
import pandas as pd
import subprocess as sp
import locate_funcs as lf
import locate_files as fi
from ew import findEW
import sys


def sigcells(linesfile, ewcut, codeLoc, flog, wave, testing=0):
    '''
    Determines which cells are the significiant contributers
    to the EW measurement. Cells are removed until the EW of 
    the absorption feature (with zero SNR) is different from the 
    EW of the absorption (0 SNR) from the full, original list of 
    cells by a percentage deteremined by ewcut
    '''

    print('In sigcells: ',wave)
    singleCellCount = 0       # Counts number of LOS dominated by a single cell
    
    # Get the info from the filename
    galID  = linesfile.split('.')[0]
    ion    = linesfile.split('.')[1]
    losnum = (linesfile.split('.')[2]).split('los')[1]

    flog.write('\n{0:s}\n'.format(losnum))

    specCommand = codeLoc+'/funcs/mkspec/specsynth los_single.list Mockspec_0SNR.runpars'
    specFileBase = '{0:s}.{1:s}.los{2:s}.{3:s}.spec'
    specfile = specFileBase.format(galID, ion, losnum, wave)

    # Get the velocity limits, EW, and vbar from sysabs file
    negVelLimit,posVelLimit,ewSysabs,vbar = lf.vel_limits(linesfile)
    ewVelcut = ewSysabs

    # Read in the linesfile
    header = 'z logN b id vpec'.split()
    fname = linesfile+'.velcut'
    with open(fname,'r') as f:
        redshift = float(f.readline().strip())
    linesdf = pd.read_csv(fname,sep='\t',names=header,skiprows=1)

    numcells = len(linesdf)

    if numcells>1:

        # Generate weights to determine the relative importance of each cell
        # based on the column density and vpec
        linesdf['vNorm'] = (linesdf['vpec']-vbar) / linesdf['vpec'].std()
        linesdf['weight'] = linesdf['logN'] + np.abs(linesdf['vNorm'])
    
        # Sort the cells by decreasing column density
        linesdf.sort_values(by='weight',ascending=False,inplace=True)

        # Create a los.list file containing only this LOS
        with open('los_single.list', 'w') as fLos:
            datfile = linesfile.replace('lines', 'dat') + '\n'
            fLos.write(datfile)
        
        # Rename the velcut .lines to remove velcut from name, 
        # so it will be used by specsynth
        command = 'cp {0:s}.velcut {0:s}'.format(linesfile)
        sp.call(command, shell=True)

        # Run specsynth on the velcut lines list
        sp.call(specCommand, shell=True)

        # Get the EW of this noise-less spectra
        specFileBase = '{0:s}.{1:s}.los{2:s}.{3:s}.spec'
        specfile = specFileBase.format(galID, ion, losnum, wave)

        # Copy the initial quiet spectra
        command = 'cp {0:s} {0:s}.velcutclean'.format(specfile)
        sp.call(command, shell=True)
        
        # Read in the spectra data
        specdata = np.loadtxt(specfile)
        wavelength = specdata[:,0]
        velocity = specdata[:,1]
        flux = specdata[:,2]
        
        # Get the EW of the full lines file with S/N = 0
        ew = findEW(wavelength, velocity, flux, negVelLimit, posVelLimit)
        ewdiff = abs( (ewSysabs - ew) / ewSysabs )
        ew_velcut_lines = ew
        ewVelcut = ew

        # See how many cells are in the velcut file
        with open(linesfile+'.velcut') as fcut:
            for i,l in enumerate(fcut):
                pass
        numcells = i   # One less than the actual number of lines due to the
                       # absorption redshift in the first line
        
    # If there is only one cell in the file, stop. Clearly that cell
    # is responsible for all the absorption
    # Write the significant lines to file
    s = '{0:>8.7f}\t{1:>8f}\t{2:>8f}\t{3:>8d}\t{4:>8f}\n'
    finalLinesFile = linesfile.replace('.lines','.{0:s}.lines.final'.format(wave))
    with open(finalLinesFile, 'w') as f:
        f.write('{0:.16f}\n'.format(redshift))

        if numcells == 1:
            singleCellCount += 1
            sigEnd = 1    
            try:
                f.write(s.format(linesdf['z'].iloc[0],
                        linesdf['logN'].iloc[0],
                        linesdf['b'].iloc[0],
                        int(linesdf['id'].iloc[0]),
                        linesdf['vpec'].iloc[0]))
            except TypeError:
                print('TypeError in significant.py: ')
                print(losnum, len(cell_ID), len(cell_z), numcells)
                print(type(cell_ID), type(cell_z))
                sys.exit()

        else:
            start = 0
            end = len(linesdf)
            depth = 0
            sigEnd = search(start, end, ewcut, linesdf, 
                redshift, specCommand, linesfile, specfile, negVelLimit, 
                posVelLimit, ewVelcut, flog, depth)

            for i in range(0,sigEnd):
                f.write(s.format(linesdf['z'].iloc[i],
                        linesdf['logN'].iloc[i],
                        linesdf['b'].iloc[i],
                        int(linesdf['id'].iloc[i]),
                        linesdf['vpec'].iloc[i]))
    # Return the wavelength used
    


def search(start, end, ewcut, linesdf,redshift, 
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
        cutdf = linesdf[:index]
        #cut_z = lines_z[:index]
        #cut_N = lines_N[:index]
        #cut_b = lines_b[:index]
        #cut_ID = lines_ID[:index]

        # Write the new array to the lines file
        s = '{0:>8.7f}\t{1:>8f}\t{2:>8f}\t{3:>8d}\t{4:>8f}\n'
        with open(linesfile, 'w') as f:            
            f.write('{0:.16f}\n'.format(redshift))
            for i in range(0,len(cutdf)):
                f.write(s.format(cutdf['z'].iloc[i],
                        cutdf['logN'].iloc[i],
                        cutdf['b'].iloc[i],
                        int(cutdf['id'].iloc[i]),
                        cutdf['vpec'].iloc[i]))

        # Run specsynth
        sp.call(specCommand, shell=True)

        # Find the new EW of the cut down list of cells
        wavelength, velocity, flux = np.loadtxt(specfile, usecols=(0,1,2),
                                                unpack=True)
        ew = findEW(wavelength, velocity, flux, negVelLimit, posVelLimit)
        ewdiff = abs( (ewVelcut - ew) / ewVelcut)*100
        flog.write('{0:d}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\t{5:.3f}\t{6:.3f}\n'.format(
                    len(linesdf), start, end, mid, len(cutdf), ew, ewdiff))
    
        if ewdiff<ewcut:
            # Not deep enough cut
            # Only send in top half
            end = mid
            search(start, end, ewcut, linesdf, 
                    redshift, specCommand, linesfile, specfile, negVelLimit, 
                    posVelLimit, ewVelcut, flog, depth)
        else:
            # Cut is too deep
            # Need to include more of cut cells 
            start = mid
            search(start, end, ewcut, linesdf, 
                    redshift, specCommand, linesfile, specfile, negVelLimit, 
                    posVelLimit, ewVelcut, flog, depth)
    return end

