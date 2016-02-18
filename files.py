

# Functions pertaining to creating, reading files and directories

# Functions contained:
#   read_control_file
#   setup_galprops
#   get_transition_info
#   setup_idcells
#   setup_mockspec
#   setup_ion_dir

import sys
import numpy as np
import os
import subprocess as sp

def read_control_file():

    """
    Read in the control file, named mockspec.config
    Returns everything in the files
    """

    filename = 'mockspec.config'
    try:
        f = open(filename)
    except IOError:
        print 'Missing control file: {0:s}'.format(filename)
        print 'Quitting...'
        sys.exit()

    # Read past header
    f.readline()
    f.readline()

    galID = f.readline().split()[0]
    expn = f.readline().split()[0]
    nlos = int(f.readline().split()[0])
    maximpact = float(f.readline().split()[0])
    incline = float(f.readline().split()[0])
    ewcut = float(f.readline().split()[0])
    snr = float(f.readline().split()[0])
    ncores = int(f.readline().split()[0])
    rootLoc = f.readline().split()[0]
    sigcellsCut = float(f.readline().split()[0])
    # Now at flags for running the various subfunctions
    f.readline()
    runRates = int(f.readline().split()[0])
    runGenLOS = int(f.readline().split()[0])
    runCellfinder = int(f.readline().split()[0])
    runIdcells = int(f.readline().split()[0])
    runLos7 = int(f.readline().split()[0])
    runSpecsynth = int(f.readline().split()[0])
    runSysanal = int(f.readline().split()[0])
    runCullabs = int(f.readline().split()[0])
    runLocateCells = int(f.readline().split()[0])
    # Now at ion section
    # Loop over rest of file
    ions = []
    xh = []
    instruments = []
    f.readline()
    f.readline()
    for line in f:
        if line.startswith('#')==False:
            l = line.split()
            ions.append(l[0])
            xh.append(l[1])
            instruments.append(l[2])

    f.close()

    props = (galID, expn, nlos, maximpact, incline, ewcut, snr, ncores, rootLoc, sigcellsCut)
    flags = (runRates, runGenLOS, runCellfinder, runIdcells, runLos7, 
             runSpecsynth, runSysanal, runCullabs, runLocateCells)
    return props, flags, ions, xh, instruments


def setup_galprops(galID, expn, requiredLoc, summaryLoc):

    from subprocess import call

    if not os.path.exists('gal_props.dat'):
        command = 'cp ' +requiredLoc+ '/gal_props.dat .'
        call(command, shell=True)

    call('cp gal_props.dat gal_props.dat.tmp', shell=True)
    f = open('gal_props.dat.tmp')
    gpf = open('gal_props.dat', 'w')

    for i in range(0,4):
        line = f.readline()
        gpf.write(line)

    gasfile = galID+'_GZa'+expn+'.txt'
    line = f.readline()
    gpf.write(line.replace('MW9_GZ932.a1.001.HI.txt', gasfile))
    
    line = f.readline()
    gpf.write(line.replace('MW9', galID))

    line = f.readline()
    gpf.write(line.replace('1.001', expn))

    gpf.write(f.readline())
    gpf.write(f.readline())
    
    gpf.write(summaryLoc)
    f.close()
    gpf.close()

    call('rm gal_props.dat.tmp', shell=True)



def get_transition_info(ion, codeLoc):

    # Returns the atomic number and existion level 
    # of this ion
    # Open the Mockspec.transistions file
    f = open(codeLoc+'/data/Mockspec.transitions')

    for line in f:
        fion = line.split()[4]
        if fion==ion:
            element = line.split()[1]
            k = line.split()[2]
            j = line.split()[3]
            f.close()
            return element, k, j




def setup_idcells(gasfile, ion_list):
    
    from subprocess import call

    if not os.path.exists('idcells_props.dat'):
        command = 'cp /home/matrix3/jrvander/requiredfiles/idcells_props.dat .'
        call(command, shell=True)

    call('cp idcells_props.dat idcells_props.tmp', shell=True)
    tmpf = open('idcells_props.tmp')
    datf = open('idcells_props.dat', 'w')

    for i in range(0,3):
        line = tmpf.readline()
        datf.write(line)

    line = tmpf.readline()
    datf.write(line.replace('dwarf9o_GZa1.001.txt', gasfile))
    
    line = tmpf.readline()
    datf.write(line)

    line = tmpf.readline()
    ion_string = ''
    for ion in ion_list:
        ion_string = ion_string+ion+' '
    datf.write(ion_string + '\t\t# Ions to use')
    tmpf.close()
    datf.close()

    call('rm idcells_props.tmp', shell=True)






def setup_mockspec(ion_list, instr_list, ewcut, snr, xh_list, requiredLoc):

    from subprocess import call

    # Need Mockspec.runpars, Mockspec.instruments, Mockspec.transitions
    basecommand = 'cp '+requiredLoc
    if not os.path.exists('Mockspec.runpars'):
        command = basecommand + 'Mockspec.runpars .'
        call(command, shell=True)
    if not os.path.exists('Mockspec.instruments'):
        command = basecommand + 'Mockspec.instruments .'
        call(command, shell=True)
    if not os.path.exists('Mockspec.transtions'):
        command = basecommand + 'Mockspec.transitions .'
        call(command, shell=True)
        
    # Alter Mockspec.runpars
    # The others should never be altered
    call('cp Mockspec.runpars Mockspec.tmp', shell=True)
    tmpf = open('Mockspec.tmp')
    datf = open('Mockspec.runpars', 'w')

    line = tmpf.readline()
    datf.write(line)

    slev = '5.0'  # Significance level of detection
    
    for ion, inst, xh in zip(ion_list, instr_list, xh_list):

        element, Z, stage = get_transition_info(ion, requiredLoc)
        if element =='hydrogen':
            vmax = '10000.'
        else:
            vmax = '1000.'
        element = '\''+element+'\''
        line = ('{0:<11s} {1:<5s} {2:<6s} {3:<7s} '
                '{4:<6f} {5:<5f} {6:<6s} {7:<s}\n'.format(element, 
                stage, xh, inst, ewcut, snr, vmax, slev))

        datf.write(line)
        
    call('rm Mockspec.tmp', shell=True)
    tmpf.close()
    datf.close()


def get_transition_info(ion, requiredLoc):

    """
    Returns the atomic number and exiciation level
    of this ion as defined in Mockspec.transitions
    """
    filename = requiredLoc+'/Mockspec.transitions'
    f = open(filename)
    for line in f:
        fion = line.split()[4]
        if fion==ion:
            element = line.split()[1]
            k = line.split()[2]
            j = line.split()[3]
            f.close()
            return element, k, j
    


def setup_ion_dir(ion, galID, expn, codeLoc):
    
    """ 
    Creates a directory for this ion if one does not exit
    Copies the all required files into this directory
    Returns the path to the ion directory
    """
    cwd = os.getcwd()
    ionloc = cwd + '/' + ion
    if not os.path.exists(ionloc):
        command = 'mkdir '+ion
        try:
            sp.check_call(command, shell=True)
        except:
            print 'Could not complete {0:s}'.format(command)

    # Create the los.list file
    command = 'ls *'+ion+'*los*dat > qso.list && mv qso.list ./'+ion+'/'
    try:
        sp.check_call(command, shell=True)
    except:
        print 'Could not complete {0:s}'.format(command)
    
    # Copy the cell files, ion boxes and the lines files into the ion directory
    ionbox =  galID+'_GZa'+expn+'.'+ion+'.txt'
    command = 'cp '+ionbox+' ./'+ion+'/'
    try:
        sp.call(command, shell=True)
    except: 
        print 'Could not complete {0:s}'.format(command)

    command = 'mv *'+ion+'.los*.dat ./'+ion+'/'
    try:
        sp.check_call(command, shell=True)
    except:
        print 'Could not complete {0:s}'.format(command)
        
    command = 'cp lines* ./'+ion+'/'
    try:
        sp.check_call(command, shell=True)
    except:
        print 'Could not complete {0:s}'.format(command)
   

    # Copy the Mockspec files from the parent directory to here
    command = 'cp '+codeLoc+'/controls/Mockspec* ./'+ion+'/'
    try:
        sp.check_call(command, shell=True)
    except:
        print 'Could not complete {0:s}'.format(command)

    # Check that nothing exists in the parent directory
    command = 'ls *{0:s}*los*dat | wc -l'.format(ion)
    numDatFiles = int(sp.check_output(command, shell=True).strip())
    if numDatFiles!=0:
        print 'ERROR in setup_ion_dir in files.py'
        print 'Problem moving .dat file for {0:s}'.format(ion)
        print 'Exitting....'
        sys.exit()
    return ionloc


def setup_inclination_dir(incline, ions, runRates, galID, expn):

    """
    Sets up a directort called i<incline> in the expansion factor
    directory. Moves all rates output files into it if rates
    was run. Returns the file path to the new directory
    """
    inc = int(incline)    

    # Check to see if the inclination directory already
    # exits
    if not os.path.exists('./i{0:d}'.format(inc)):
        command = 'mkdir i{0:d}'.format(inc)
        sp.call(command, shell=True)

    cwd = os.getcwd()
    incLoc = '{0:s}/i{1:d}/'.format(cwd, inc)

    # Move the ion files output by rates into this
    # new directory
    if runRates==1:
        for ion in ions:
            boxName = '{0:s}_GZa{1:s}.{2:s}.txt'.format(galID, expn, ion)
            command = 'mv {0:s} ./i{1:d}/'.format(boxName, inc)
            try:
                sp.check_call(command, shell=True)
            except CalledProcessError:
                print 'Could not complete {0:s} in setup_inclination_dir'.format(command)
                sys.exit()
    else:
        # Ensure the ion boxes are in the directory
        for ion in ions:
            boxName = '{0:s}_GZa{1:s}.{2:s}.txt'.format(galID, expn, ion)
            boxPath = incLoc+boxName

            # Check if the ion box is in the inclination directory already
            if not os.path.isfile(boxPath):
                # Ion boxes are not in the inclination directory
                # Check if they are in the expansion directory
                newPath = '{0:s}/{1:s}'.format(cwd, boxName)
                if os.path.isfile(newPath):
                    # Move the box to the inclination directory
                    command = 'mv {0:s} ./i{1:d}/'.format(boxName, inc)
                    try:
                        sp.check_call(command, shell=True)
                    except CalledProcessError:
                        print 'Could not complete {0:s} in setup_inclination_dir'.format(command)
                else:
                    # Ion box is missing
                    print 'Cannot find ion box in setup_inclination_dir'
                    print 'Exiting...'
                    sys.exit()

    # Copy the rest of the control files into the directory
    boxname = '{0:s}_GZa{1:s}.txt'.format(galID, expn)
    filenames = ['mockspec.config', 'gal_props.dat', 'galaxy.props', boxname]
    for fn in filenames:
        command = 'cp {0:s} ./i{1:d}/'.format(fn, inc)
        try:
            sp.check_call(command, shell=True)
        except CalledProcessError:
            print 'Cannot find {0:s} in setup_inclination_dir'.format(fn)
            sys.exit()


    # Return the path to the new directory
    return incLoc


def setup_galaxy_props(summaryLoc, galID, expn, inc):
        
    """
    Creates a file called galaxy.props in the root directory.
    Contains relevent galaxy information from the <galID>.dat file
    stored in the codeLoc's summary directory
    """
    
    cwd = os.getcwd()
    # Read in properties from the summary file
    summaryFile = '{0:s}.dat'.format(galID)
    try:
        f = open(summaryLoc+summaryFile)
    except IOError:
        try:
            f = open('../output/rotmat_a{0:s}.txt'.format(expn))
            except IOError:
                print 'Could not open summary file:\n\t{0:s}{1:s}'.format(summaryLoc, summarFile)
                print 'Coult not find rotmat file'
                sys.exit()
    f.readline()
    f.readline()
    found = 0
    for line in f:
        l = line.split()
        a = l[0]
        if a==expn:
            found = 1
            redshift = float(l[1])
            mvir = float(l[2])
            rvir = float(l[3])
    if found==0:
        print 'Could not find {0:s} in {1:s}{2:s}'.format(expn, summaryLoc, summaryFile)
        sys.exit()
    f.close()

    f = open('galaxy.props', 'w')

    s = 'galID          {0:s}\n'.format(galID)
    f.write(s)
    s = 'Expn           {0:s}\n'.format(expn)
    f.write(s)
    s = 'Redshift       {0:.3f}\n'.format(redshift)
    f.write(s)
    s = 'Mvir           {0:.4e}\n'.format(mvir)
    f.write(s)
    s = 'Rvir           {0:.4f}\n'.format(rvir)
    f.write(s)
    s = 'Inclination    {0:.1f}\n'.format(inc)
    f.write(s)

    f.close()


    


def rename(galID, expn, ion, incline, runLocateCells, runCullabs):

    '''
    Renames files output by the pipeline to incline the 
    inclination angle
    '''
    inc = int(incline)
    if runCullabs==1:
        # Need to rename ALL.sysabs file
        allName = '{0:s}.{1:s}.a{2:s}.ALL.sysabs'.format(galID, ion, expn)
        newName = allName.replace('ALL','i{0:d}.ALL'.format(inc))
        command = 'mv {0:s} {1:s}'.format(allName, newName)
        try:
            sp.check_call(command, shell=True)
        except:
            print 'Error in rename running \n\t{0:s}'.format(command)
            print 'Exitting...'
            sys.exit()

    if runLocateCells==1:
        # Need to rename the abs_cells file
        absName = '{0:s}.{1:s}.{2:s}.abs_cells.dat'
        newName = absName.replace('abs','i{0:d}.abs'.format(inc))
        command = 'mv {0:s} {1:s}'.format(absName, newName)
        try:
            sp.check_call(command, shell=True)
        except:
            print 'Error in rename running \n\t{0:s}'.format(command)
            print 'Exitting...'
            sys.exit()

       








