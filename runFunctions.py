
# Functions to run external codes as part of the pipeline

import sys
from subprocess import check_call

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
        check_call(command, shell=True)
    except:    
        print '\n\nCould not run los7 with :\n\t{0:s}'.format(command)
        print 'Exiting...'
        sys.exit()
     


def specsynth(codeLoc):

    """
    Runs specsynth, which generates the mock spectra
    Generates two .spec files for each LOS
    Takes in the local location of the code and 
    returns nothing.
    """
    
    funcLoc = '/funcs/mkspec/specsynth'
    command = codeLoc + funcLoc

    
    try:
        check_call(command, shell=True)
    except:    
        print '\n\nCould not run specsynth with :\n\t{0:s}'.format(command)
        print 'Exiting...'
        sys.exit()
     

def sysanal(codeLoc):

    """
    Runs sysanal, which calculates the EW and AOD column density for
    the mock spectra.
    Generates a .sysabs file for each LOS
    Takes in the local location of the code, returns nothing.
    """
    
    funcLoc = '/funcs/anal/sysanal'
    command = codeLoc + funcLoc

    
    try:
        check_call(command, shell=True)
    except:    
        print '\n\nCould not run sysanal with :\n\t{0:s}'.format(command)
        print 'Exiting...'
        sys.exit()
     
    


def cullabs(codeLoc):

    """
    Runs cullabs, which summarizes the results of sysanal into 
    an ALL.sysabs file
    Takes in the local location of the code and 
    returns nothing.
    """
    
    funcLoc = '/funcs/mkALL/cullabs'
    command = codeLoc + funcLoc

    
    try:
        check_call(command, shell=True)
    except:    
        print '\n\nCould not run cullabs with :\n\t{0:s}'.format(command)
        print 'Exiting...'
        sys.exit()
     

