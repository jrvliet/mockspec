

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

    # Reads in the control file
    # Filename: mockspec.config
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
    requiredLoc = f.readline().split()[0] 
    # Now at ion section
    # Loop over rest of file
    ions = []
    xh = []
    instruments = []
    f.readline()
    f.readline()
    for line in f:
        l = line.split()
        ions.append(l[0])
        xh.append(l[1])
        instruments.append(l[2])

    f.close()

    props = (galID, expn, nlos, maximpact, incline, ewcut, snr, ncores, rootLoc, requiredLoc)

    return props, ions, xh, instruments


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






def setup_mockspec(ion_list, instr_list, ewcut, snr, sh_list):

    from subprocess import call

    # Need Mockspec.runpars, Mockspec.instruments, Mockspec.transitions

    if not os.path.exists('Mockspec.runpars'):
        command = 'cp /home/matrix3/jrvander/requiredfiles/Mockspec.runpars .'
        call(command, shell=True)
    if not os.path.exists('Mockspec.instruments'):
        command = 'cp /home/matrix3/jrvander/requiredfiles/Mockspec.instruments .'
        call(command, shell=True)
    if not os.path.exists('Mockspec.transtions'):
        command = 'cp /home/matrix3/jrvander/requiredfiles/Mockspec.transitions .'
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

        element, Z, stage = getTransitionInfo(ion)
        if element =='hydrogen':
            vmax = '5000.'
        else:
            vmax = '1000.'
        element = '\''+element+'\''
        line = '{0:<11s} {1:<5s} {2:<6s} {3:<7s} {4:<6s} {5:<5s} {6:<6s} {7:<s}\n'.format(
element, stage, xh, inst, ewcut, snr, vmax, slev)
        datf.write(line)
        
    call('rm Mockspec.tmp', shell=True)
    tmpf.close()
    datf.close()




def setup_ion_dir(ion, galID, expn):
    
    """ 
    Creates a directory for this ion if one does not exit
    Copies the all required files into this directory
    Returns the path to the ion directory
    """
    cwd = os.getcwd()
    ionloc = cwd + '/' + ion
    if not os.path.exists(loc):
        command = 'mkdir '+ion
        sp.call(command, shell=True)

    # Copy the cell files, ion boxes and the lines files into the ion directory
    ionbox =  galID+'_GZa'+expn+'.'+ion+'.txt'
    command = 'cp '+ionbox+' ./'+ion+'/'
    sp.call(command, shell=True)
    command = 'mv *'+ion+'.los*.dat ./'+ion+'/'
    sp.call(command, shell=True)
    command = 'cp lines* ./'+ion+'/'
    sp.call(command, shell=True)
   
     # Create the los.list file
    command = 'ls *los*dat > los.list'
    sp.call(command, shell=True)

    # Copy the Mockspec files from the parent directory to here
    command = 'cp ../Mockspec* .'
    sp.call(command, shell=True)

    return ionloc




