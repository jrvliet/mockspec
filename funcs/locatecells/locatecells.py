# 
# Description:
#  Determines which cells are the significant contributers 
#  to the EW measurement for all LOS for a given ion/snapshot
#

from __future__ import print_function
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


def locateSigCells(run,ion,codeLoc,testing=0):

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
    boxfile = '{0:s}_GZa{1:s}.{2:s}.h5'.format(run.galID,run.expn,ion.name)
    box = pd.read_hdf(boxfile, 'data')
    if testing==1:
        print('Box read in')

    # Read in the LOS info from lines.info
    los_info = np.loadtxt('lines.info',skiprows=2)

    # Determine which transition to use
    # Only need on, use the transition with the weaker oscillator strength since
    # this will retain mroe cells than the other transitions
    waves = fi.transition_name(ion.name,testing)

    # Open the output file
    outfile = '{0:s}.{1:s}.{2:s}.i{3:d}.abs_cells.h5'.format(run.galID,
                        run.expn,ion.name,int(run.incline))
    header = ['wave','LOS','D','cellID','redshift','logN','dobbler_b',
                'x', 'y', 'z', 'vx', 'vy', 'vz',
                'r', 'nH', 'temperature', 'cell_size', 'SNII', 'SNIa', 
                'alpha_Zmet', 'ion_density', 'fullLogNstart', 'fullLogNend']

    # Open the summary file
    sumfile = '{0:s}_a{1:s}_{2:s}_i{3:d}_absCellsSummary.h5'.format(run.galID,
                        run.expn,ion.name,int(run.incline))
    sumheader = ['los', 'impact', 'phi', 'wave','startNumCells', 'startEW', 'startlogN', 
                'endNumCells', 'endEW', 'endlogN']
    summary = pd.DataFrame(columns=sumheader)

    # Write a header to the output file
    numcols = len(header)

    # Create a blank row of zeros to build the array with
    d = pd.DataFrame(columns=header)

    # Make a version of Mockspec.runpars that has zero SNR
    # Needed for sigcells
    lf.quiet_mockspec()

    # Get a list of LOS that have a sysabs file associated with it
    sysabs_losnum = glob.glob('*los*sysabs')
    sysabs_losnum.sort()
    
    if testing==1:
        print('Sysabs files aggregated')
        print('Number of sysabs files: ', len(sysabs_losnum))

    flog = open('sig_cells.log', 'w')
    falselog = open('falseDetections.log','w')
    falselog.write('los\tEWsysabs\tv-\tv+\tInitNumCells\n')

    # Loop over lines of sight
    for i in range(0,len(sysabs_losnum)):

        losnum = sysabs_losnum[i].split('.')[2].split('los')[1]
        num = int(losnum)

        linesfile = '{0:s}.{1:s}.los{2:s}.lines'.format(run.galID,ion.name,losnum)

        # Get the column density of the LOS from the lines file
        logNinitial = lf.linesLogN(linesfile)

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
        phi = los_info[num-1,2]
        
        # Loop over the different transitions
        for wave in waves:
            losSummary = pd.Series(index=sumheader)
            losSummary['wave'] = wave
            losSummary['los'] = num
            losSummary['startNumCells'] = numCells
            losSummary['impact'] = imp
            losSummary['phi'] = phi
            losSummary['startlogN'] = lf.linesLogN(linesfile)

            # Copy the original lines file
            command = 'cp {0:s} {0:s}.tmp'.format(linesfile)
            sp.call(command, shell=True)

            # Perform the velocity cut
            lf.velcut(linesfile, testing=testing)

            # Find the significant cells
            endCut,startEW,endEW = sg.sigcells(linesfile,run.sigcellsCut,codeLoc,
                                                flog,falselog,wave,testing=testing)

            # Get the properties of the cells
            # Open the lines.final file
            finalLinesFile = linesfile.replace('.lines',
                                '.{0:s}.lines.final'.format(wave))
            logNfinal = lf.linesLogN(finalLinesFile)
            losSummary['endlogN'] = logNfinal
            losSummary['startEW'] = startEW
            losSummary['endEW'] = endEW

            final_file = open(finalLinesFile)
            final_file.readline()

            endNumCells = 0
    
            # Loop over the significant cells
            for line in final_file:
                endNumCells += 1
                l = line.split()
                redshift = float(l[0])
                column = float(l[1])
                doppler = float(l[2])
                cellID = int(float(l[3]))
            
                # Get the cell's properties from the boxfile
                index = cellID-1
                x = box['x'].iloc[index]
                y = box['y'].iloc[index]
                z = box['z'].iloc[index]
            
                # Calculate the galactocentric distance
                r = np.sqrt(x*x + y*y + z*z)
                
                cell = pd.Series(index=header)
                cell['wave'] = wave
                cell['LOS'] = num
                cell['D'] = imp
                cell['cellID'] = cellID
                cell['redshift'] = redshift
                cell['logN'] = column 
                cell['dobbler_b'] = doppler
                cell['x'] = x
                cell['y'] = y 
                cell['z'] = z 
                cell['vx'] = box['vx'].iloc[index]
                cell['vy'] = box['vy'].iloc[index]
                cell['vz'] = box['vz'].iloc[index]
                cell['r'] = r
                cell['nH'] = np.log10(box['nH'].iloc[index])
                cell['temperature'] = np.log10(box['temperature'].iloc[index])
                cell['cell_size'] = box['cell_size'].iloc[index]
                cell['SNII'] = np.log10(box['SNII'].iloc[index])
                cell['SNIa'] = np.log10(box['SNIa'].iloc[index])
                cell['alpha_Zmet'] = box['alpha_Zmet'].iloc[index]
                cell['ion_density'] = np.log10(box['nIon'].iloc[index])
                cell['fullLogNstart'] = logNinitial
                cell['fullLogNend'] = logNfinal
           
                # Append to the main array
                d = d.append(cell,ignore_index=True)
                
            losSummary['endNumCells'] = endNumCells
            summary = summary.append(losSummary,ignore_index=True)

            # Rename the original .lines file back to its full name
            command = 'cp {0:s}.tmp {0:s}'.format(linesfile)
            sp.call(command, shell=True)
          
    # Write the outfile
    d.to_hdf(outfile,'data',mode='w')
    summary.to_hdf(sumfile,'data',mode='w')

    print('For {0:s}, {1:d} LOS are dominated by one cell'.format(ion.name,
            singleCount))

    flog.close()
    falselog.close()
