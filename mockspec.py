#
#  mockspec.py
# 
#  Full pipeline for performing synthetic QSO observations
#  of ART simulations of galaxy CGM

# General libraries
import numpy as np
import tables as tb
import os
import pandas
import subprocess as sp

# Code specific libraries
from files import *
from genLOS import *
from ratesControl import *



pathname = os.path.dirname(sys.argv[0])
codeLoc = os.path.abspath(pathname)
summaryLoc = codeLoc+'/summaries/'


#  Read in the control file
print 'Reading in control file...'
props, ions, xh, instruments = readControlFile()
galID, expn, nlos, maximpact, incline, ewcut, snr, ncores, rootLoc, requiredLoc = props

# Genearte the name of the gasfile
gasfile = galID+'_GZa'+expn+'.txt'
 
print 'Gasfile:      ', gasfile
print 'Ions:         ', ions
print 'NLOS:         ', nlos
print 'Max impact:   ', maximpact
print 'Incline:      ', incline
print 'EWCut:        ', ewcut
print 'SNR:          ', snr
print 'NCores:       ', ncores
print 'Root Loc:     ', rootLoc
print 'Requried Loc: ', requiredLoc

# Generate gal_props.dat file
print 'Generating gal_profgrep ps.dat...'
setupGalprops( galID, expn, requiredLoc, summaryLoc )


##### 
#  
#  Run rates
#
#####
print '\nRates:'
print '\t Setting up rates.inp...'
setupRatesControl( gasfile, expn, ions, requiredLoc)
print '\t Setting up rates.outfiles...'
setupRatesOutputs(galID, expn, ions, codeLoc, requiredLoc) 
print '\t Running rates...'
#runRates(codeLoc)


##### 
#  
#  Generate the lines of sight
#
#####
print '\nGenerating lines of sight...'
genLines(galID, gasfile, summaryLoc, expn, incline, nlos, maximpact, ncores)

##### 
#  
#  Run lines of sight through box with cellfinder
#
#####
print '\nRunning LOS through box...'
runCellfinder(codeLoc)



# Start looping over ions
for ion in ions:


    # Setup the ion directory
    ionloc = setup_ion_dir(ion, galID, expn) 

    # Move into the ion directory
    os.chdir(ionloc)

    ##### 
    #  
    #  los7
    #  Determine path length through each cell and create the
    #  .losdata and .lines files
    #
    #####
    print '\nDetermining cell path length and applying rough cut'
    los7(codeLoc)


    ##### 
    #  
    #  specsynth
    #  Generate the synthetic spectra
    #
    #####
    print '\nGenerating spectra'
    specsynth(codeLoc)


    ##### 
    #  
    #  sysanal
    #  Analyze the spectra
    #
    #####
    print '\nAnalyzing spectra...'
    sysanal(codeLoc)



    ##### 
    #  
    #  cullabs
    #  Create sysabs file
    #
    #####
    print '\nCreating sysabs'
    sysabs(codeLoc)

    

    ##### 
    #  
    #  locatecells
    #  Identify the cells that are significant contributers to the absorption
    #
    #####
    print '\nIdentifying significant cells'
    runLocateCells(codeLoc)

    os.chdir('..')





