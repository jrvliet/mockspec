
import numpy as np
import subprocess as sp
import locate_funcs as lf
from ew import findEW

def sigcells(linesfile, ewcut, codeLoc, testing=0):
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

    # Read in the linesfile
    cell_z, cell_N, cell_b, cell_ID = np.loadtxt(linesfile+'.velcut', skiprows=1,
                                    usecols=(0,1,2,3), unpack=True)

    if testing==1:
        print 'In sigcells, number of velcut cells read in: ', len(cell_z)

    # Get the EW from sysabs
    negVelLimit, posVelLimit, ewSysabs = lf.vel_limits(linesfile)

    #################################################################
    #                                                               #
    #     Generate a noise-less spectra for the velcut .lines file  #
    #                                                               #
    #################################################################

    # Create a los.list file containing only this LOS
    with open('los_single.list', 'w') as fLos:
        datfile = linesfile.replace('lines', 'dat') + '\n'
        fLos.write(datfile)
    
    # Rename the velcut .lines to remove velcut from name, 
    # so it will be used by specsynth
    command = 'cp '+linesfile+'.velcut '+linesfile
    sp.call(command, shell=True)

    # Run specsynth on the velcut lines list
    specsynth_command = codeLoc+'/funcs/mkspec/specsynth los_single.list Mockspec_0SNR.runpars'
    sp.call(specsynth_command, shell=True)

    # Get the EW of this noise-less spectra
    bluewave, redwave = lf.transition_name(ion, codeLoc)

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
    if numcells == 1:
        singleCellCount += 1
    
    else:
        # More than one cell in the .lines.velcut file
        # Find the ones that are significant

        # Read in .lines.velcut file
        with open(linesfile+'.velcut') as fvelcut:
            redshift = float(fvelcut.readline().strip())

        velcutdata = np.loadtxt(linesfile+'.velcut', skiprows=1)

        velcut_z = list(velcutdata[:,0])
        velcut_N = list(velcutdata[:,1])
        velcut_b = list(velcutdata[:,2])
        velcut_ID = list(velcutdata[:,3])
        
        f_log = open('sigcells.log', 'w')
        f_log.write('Numcells \t EW \t EWdiff\n')
        sFormat = '{0:d} \t \t{1:0.3f} \t {2:0.3f}\n'
        f_log.write(sFormat.format(len(velcut_z), ew, ewdiff))

        ewdiff = 0.5*ewcut
        loopcount = 0
        while ewdiff < ewcut:
            loopcount += 1

            # Get the number of cells
            numcells = len(velcut_z)
            
            # Cut this in half to get the midpoint
            cutIndex = int(numcells/2)

            # Only work with the top half of the list
            cut_z = velcut_z[:cutIndex]
            cut_N = velcut_N[:cutIndex]
            cut_b = velcut_b[:cutIndex]
            cut_ID = velcut_ID[:cutIndex]

            # Write the new array to the lines file
            s = '{0:>8.7f}\t{1:>8f}\t{2:>8f}\t{3:>8d}\n'
            with open(linesfile, 'w') as f:            
                f.write('{0:.16f}\n'.format(redshift)
                for i in range(0,len(cut_z)):
                    f_newlines.write(s.format(cut_z[i], cut_N[i],
                                              cut_b[i], int(cut_ID[i])))
            # Run specsynth
            sp.call(specsynth_command, shell=True)

            # Find the new EW of the cut down list of cells
            wavelength, velocity, flux = np.loadtxt(specfile, usecols=(0,1,2),
                                                    unpack=True)
            ew = findEW(wavelength, velocity, flux, negVelLimit, posVelLimit)
            ewdiff = abs( (ew_velcut_lines - ew) / ew_velcut_lines)*100
            
            # Determine if this cut is too far or not far enough
            if ewdiff > ewcut:
                # Too far
                # New midpoint is halfway between the current midpoint and
                # the end of the list
                newmid = (midpoint + len(velcut_z)) / 2.0
            else:
                # Not far enough
                # New midpoint is halfway between the current 
                newmid = midpoint / 2.0

            newcut_z = velcut_z[:,newmid]
            newcut_N = velcut_N[:,newmid]
            newcut_b = velcut_b[:,newmid]
            newcut_ID = velcut_ID[:,newmid]



            # Check that there are more than one cell left
            if len(velcut_z)>1:

                # Find the cell with the lowest column denstiy
                index = velcut_N.index(min(velcut_N))
                
                # Delete this index
                del velcut_z[index]
                del velcut_N[index]
                del velcut_b[index]
                del velcut_ID[index]

                # Write the new .lines values to file
                f_newlines = open(linesfile, 'w')
                f_newlines.write('{0:.16f}\n'.format(redshift))
                s = '{0:>8.7f}\t{1:>8f}\t{2:>8f}\t{3:>8d}\n'
                for i in range(0,len(velcut_z)):
                    f_newlines.write(s.format(velcut_z[i], velcut_N[i],
                                              velcut_b[i], int(velcut_ID[i])))
                f_newlines.close()

                # Run specsynth again
                sp.call(specsynth_command, shell=True)

                # Find the new EW 
                specdata = np.loadtxt(specfile)
                wavelength = specdata[:,0]
                velocity = specdata[:,1]
                flux = specdata[:,2]
                
                ew = findEW(wavelength, velocity, flux, negVelLimit, posVelLimit)
                ewdiff = abs( (ew_velcut_lines - ew) / ew_velcut_lines)*100
                f_log.write('{0:d} \t \t{1:0.3f} \t {2:0.3f}\n'.format(
                            len(velcut_z), ew, ewdiff))
            else:
                singleCellCount += 1
                ewdiff = ewcut
            
        f_log.close()

    # Copy the lines file for protection
    command = 'cp '+linesfile+' '+linesfile+'.final'
    sp.call(command, shell=True)

    return singleCellCount


