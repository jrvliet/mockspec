
# Functions to run external codes as part of the pipeline

import os
import sys
import subprocess as sp
#from subprocess import check_call

#def genLOS(codeLoc, galID, summaryLoc, expn, inc, nLOS, maximpact, ncores):
def genLOS(codeLoc, summaryLoc, run):
    
    """ 
    Runs genLOS, which generates random LOS for the simulation.
    Handels all rotation of the simulation to achieve the desired 
    inclination. Creates lines.dat and lines.info
    """
    funcLoc = '/funcs/generateLines/genLOS '
    args = '{0:s} {1:s} {2:s} {3:f} {4:d} {5:f} {6:d}'.format(
            run.galID, summaryLoc, run.expn, run.incline, 
            run.nlos, run.maximpact, run.ncores)
    command = codeLoc + funcLoc + args
    
    try:
        sp.run(command, shell=True)
    except Exception as e:    
        print('\n\nCould not run genLOS with :\n\t{0:s}'.format(command))
        print('Exiting...')
        print(e.message, e.args)
        sys.exit()



def runCellfinder(codeLoc, numcores):
    
    """
    Runs cellfinder, which determines which cells lie along each LOS.
    Is fully paralellzied. Creates cellID files.
    """
    cwd = os.getcwd()
    command = codeLoc+'/funcs/cellfinder/cellfinder {0:d} {1:s}'.format(numcores, cwd)
    try:
        #sp.check_call(command, shell=True)
        sp.run(command, shell=True)
    except Exception as e:
        print('\n\nCould not run cellfinder with:\n\t{0:s}'.format(command))
        print('Exiting...')
        print(e.message,e.args)
        sys.exit()


def los7(codeLoc):

    """
    Runs los7, which determines the path length of the LOS
    through each cell. Also performs a rough cut to remove 
    cells that are too small to contribute significantly
    to the final profile.
    """

    funcLoc = '/funcs/mklos/los7'
    command = codeLoc + funcLoc

    try:
        sp.check_call(command, shell=True)
    except Exception as e:
        print('\n\nCould not run los7 with :\n\t{0:s}'.format(command))
        print('In ',os.getcwd())
        print('Exiting...')
        print(e,e.args)
        sys.exit()
     


def specsynth(codeLoc):

    """
    Runs specsynth, which generates the mock spectra
    Generates two .spec files for each LOS
    Takes in the local location of the code and 
    returns nothing.
    """
    
    funcLoc = '/funcs/mkspec/specsynth qso.list'
    command = codeLoc + funcLoc

    
    try:
        sp.check_call(command, shell=True)
    except Exception as e:
        print('\n\nCould not run specsynth with :\n\t{0:s}'.format(command))
        print('Exiting...')
        print(e.message,e.args)
        sys.exit()
     

def sysanal(codeLoc):

    """
    Runs sysanal, which calculates the EW and AOD column density for
    the mock spectra.
    Generates a .sysabs file for each LOS
    Takes in the local location of the code, returns nothing.
    """
    
    funcLoc = '/funcs/anal/sysanal qso.list'
    command = codeLoc + funcLoc

    
    try:
        sp.check_call(command, shell=True)
    except Exception as e:
        print('\n\nCould not run sysanal with :\n\t{0:s}'.format(command))
        print('Exiting...')
        print(e.message,e.args)
        sys.exit()
     
    


def cullabs(codeLoc):

    """
    Runs cullabs, which summarizes the results of sysanal into 
    an ALL.sysabs file
    Takes in the local location of the code and 
    returns nothing.
    """
    
    funcLoc = '/funcs/mkALLabs/cullabs qso.list'
    command = codeLoc + funcLoc

    
    try:
        sp.check_call(command, shell=True)
    except Exception as e:
        print('\n\nCould not run cullabs with :\n\t{0:s}'.format(command))
        print('Exiting...')
        print(e.message,e.args)
        sys.exit()
     

