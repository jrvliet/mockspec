#
#  mockspec.py
# 
#  Full pipeline for performing synthetic QSO observations
#  of ART simulations of galaxy CGM

import numpy as np
import tables as tb
from mockspec_funcs import *
from genLOS import *
import subprocess as sp
from ratesControl import *
import os

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
print 'Root Loc:     ',rootLoc
print 'Requried Loc: ', requiredLoc

# Generate gal_props.dat file
print 'Generating gal_props.dat...'
setupGalprops( galID, expn, requiredLoc )


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
runRates(codeLoc)


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
command = './cellfinder'



