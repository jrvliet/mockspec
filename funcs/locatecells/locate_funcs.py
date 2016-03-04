
import numpy as np
import sys
import subprocess as sp
from ew import findEW 

def vel_limits(linesfile):
    '''
    Reads in the sysabs file corresponding to the linesfile passed in
    Returns the velocity limits and ew of the absorption 
    '''
    
    # Open the sysabs file
    sysabsfile  =  linesfile.replace('lines', 'sysabs')
    f_sysabs = open(sysabsfile)
    f_sysabs.readline()
    line = f_sysabs.readline()
    neg = float(line.split()[1])
    pos = float(line.split()[2])
    ew = float(line.split()[3])
    f_sysabs.close()

    return neg, pos, ew


def quiet_mockspec():
    '''
    Makes a version of Mockspec.runpars that has zero SNR
    '''

    # Create a Mockspec.runpars file with SNR set to zero
    f_runpars_old = open('Mockspec.runpars')
    f_runpars_new = open('Mockspec_0SNR.runpars', 'w')
    l = f_runpars_old.readline().split()
    s = '{0:s}\t\t{1:s}\t\t{1:s}\t\t{1:s}\t\t{1:s}\t\t{1:s}\t\t{1:s}\t\t{1:s}\n'
    f_runpars_new.write(s.format(l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7]))
    for line in f_runpars_old:
        l = line.split()
        l[5] = '0.'
        f_runpars_new.write(s.format(l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7]))
    f_runpars_new.close()
    f_runpars_old.close()


def velcut(linesfile, testing=0):

    """
    Cuts the cells out of the .lines file that do not fall withing
    the velocity window as found in the sysabs file
    """
    # Define constants
    c = 3.0e5   # Speed of light in km/s

    cell_z, cell_N, cell_b, cell_ID = np.loadtxt(linesfile, skiprows=1, 
                                    usecols=(0,1,2,3), unpack=True)

    # Read in the redshift of the absorption
    # This is the first line of the .lines file
    with open(linesfile, 'r') as f:
        redshift = float(f.readline().strip())

    # Get the velcoity limits of the absorption
    negVelLim, posVelLim, ewSysabs = vel_limits(linesfile)

    if testing==1:
        print '\t\tBefore velcut, number of cells: ', len(cell_z)
        print '\t\tFrom sysabs:'
        print '\t\t\tNeg Vel Limt: {0:f}'.format(neg_vel_limit)
        print '\t\t\tPos_vel_limi: {0:f}'.format(pos_vel_limit)
        print '\t\t\tEW:           {0:f}'.format(EW_sysabs)

    # New .lines file
    newlinesfile = linesfile+'.velcut'
    f_newlines = open(newlinesfile, 'w')
    
    # Write the first line
    f_newlines.write('{0:.16f}\n'.format(redshift))
    
    # Loop through the original .lines file
    velcutCount = 0
    s = '{0:>8.7f}\t{1:>8f}\t{2:>8f}\t{3:>8d}\n'
    for i in range(0,len(cell_z)):

        # Calcuate the peculiar velocity of the cell
       vpec = c*( (cell_z[i]-redshift) / (1.0+redshift) )

        # If the cell is inside the velocity range, write to file
        if vpec>neg_vel_limit and vpec<pos_vel_limit:
            f_newLines.write(s.format(cell_z[i], cell_N[i], cell_b[i], cell_ID[i])
            velcutCount += 1

    f_newlines.close()
    
    if testing==1:
        print '\t\tAfter velcut, number of cells: ', velcutCount

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

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
    negVelLimit, posVelLimit, ewSysabs = vel_limits(linesfile)

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
    if ion == 'HI':
        bluewave = 'Lya'
        redwave = ion+'1026'        
    elif ion == 'MgII':
        bluewave = ion+'2796'
        redwave = ion+'2803'
    elif ion == 'CIV':
        bluewave = ion+'1548'
        redwave = ion+'1551'
    elif ion == 'OVI':
        bluewave = ion+'1032'
        redwave = ion+'1038'
    else:
        print 'Unkown ion'
        sys.exit()

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
        for i,l in enumerate(f):
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
        velcutdata = np.loadtxt(linesfile+'.velcut', skiprows=1)
        if testing==1:
            print 'Velcutdata: ', velcutdata
            shap = velcutdata.shape
            print 'Velcutdata shape: ', shap
            print 'Number of rows: ', shap[0]
            print 'Length of shap: ', len(shap)
            print 'Index 0: ', velcutdata[1]

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
                                              velcut_b[i], velcut_ID[i]))
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


