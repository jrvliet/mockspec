import numpy as np
import sys
import subprocess as sp


def vel_limits(linesfile):
    '''
    Reads in the sysabs file corresponding to the linefile passed in
    Returns the velocity limits and ew of the absorption
    '''
    
    # Open the sysabs file
    sysabsfile = linesfile.replace('lines','sysabs')
    f = open(sysabsfile, 'r')
    f.readline()
    line = f.readline()
    neg = float(line.split()[1])
    pos = float(line.split()[2])
    ew  = float(line.split()[3])
    f.close()

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
    for line in f_runpars_old
        l = line.split()
        l[5] = '0.'
        f_runpars_new.write(s.format(l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7]))
    f_runpars_new.close()
    f_runpars_old.close()


def transition_name(ion, codeLoc):
    '''
    Finds the transition name for the ion
    Found in the Mockspec.transitions file
    '''
    waves = []
    mockspecFile = '{0:s}/controls/Mockspec.transitions'.format(codeLoc)
    with open(mockspecFile, 'r') as f:
        f.readline()
        for line in f:
            l = line.split()
            flag = int(l[0])
            ionName = l[4]
            trans = l[5]
            
            if flag==1 and ionName==ion:
                waves.append(trans)

    bluewave = waves[0]
    redwave = waves[1]

    return bluewave, redwave











