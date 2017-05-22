
import subprocess as sp
import os
import sys
from mockspec_funcs import getTransitionInfo


def setup_rates_control(gasfile, expn, ion_list, requiredLoc):

#    if not os.path.exists('rates.inp'):
    command = 'cp '+requiredLoc+'rates.inp .'
    sp.call(command, shell=True)

    # Alter the rates.inp file
    numions = len(ion_list)
    defaultratesinp = open('rates.inp')
    ratesinp = open('rates.inp.tmp', 'w')
    for i in range(0,3):
        line = defaultratesinp.readline()
        ratesinp.write(line)

    # Replace the gas file name
    line = defaultratesinp.readline().replace('MW9_GZ932.a1.001.txt', gasfile)
    ratesinp.write(line)

    # Replace the expansion parameter
    line = defaultratesinp.readline().replace('1.001', expn)
    ratesinp.write(line)

    # Replace the number of output cubes
    line = defaultratesinp.readline().replace('4', str(numions))
    ratesinp.write(line)

    # Copy the rest of the file
    for line in defaultratesinp:
        ratesinp.write(line)

    defaultratesinp.close()
    ratesinp.close()

    command = 'rm rates.inp'
    sp.call(command, shell=True)
    command = 'mv rates.inp.tmp rates.inp'
    sp.call(command, shell=True)

def setup_rates_data(codeLoc):
    '''
    Write a file with the location of the data files
    because fortran is a terrible language that can't handle
    the incredibly complicated task of command line arguements
    '''

    with open('ratesDataLoc.txt', 'w') as f:
       f.write(codeLoc+'\n')
 


def setup_rates_outputs(galID, expn, ion_list, codeLoc, requiredLoc):

    if not os.path.exists('rates.outfiles'):
        command = 'cp '+requiredLoc+'rates.outfiles .'
        sp.call(command, shell=True)
    
    # Alter the rates.outfiles file
    defoutfile = open('rates.outfiles')
    outfiles = open('rates.outfiles.tmp', 'w')
    for ion in ion_list:
        ionbox = galID+'_GZa'+expn+'.'+ion+'.txt'
        element, Z, excitation = getTransitionInfo(ion, codeLoc)
        line = '{0:<27s} {1:>2s} {2:>2s}\n'.format(ionbox, Z, excitation)
        outfiles.write(line)
    outfiles.close()
    defoutfile.close()
    command = 'rm rates.outfiles'
    sp.call(command, shell=True)
    command = 'mv rates.outfiles.tmp rates.outfiles'
    sp.call(command, shell=True)
   

def run_rates(codeLoc):
    
    command = '{0:s}/funcs/rates/rates {0:s}'.format(codeLoc)
    try:
        sp.check_call(command, shell=True) 
    except:
        print('\n\nCould not run rates with \n\t{0:s}'.format(command))
        print('Exiting....')
        sys.exit()







