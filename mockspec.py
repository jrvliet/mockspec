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
props, ions, xh, instruments = readControlFile()
galID, expn, nlos, maximpact, incline, ewcut, snr, ncores, rootLoc, requiredLoc = props

# Genearte the name of the gasfile
gasfile = galID+'_GZa'+expn+'.txt'
 
# Generate gal_props.dat file
setupGalprops( galID, expn, requiredLoc )


##### 
#  
#  Run rates
#
#####
#setupRatesControl( gasfile, expn, ions, requiredLoc)
#setupRatesOutputs(galID, expn, ions, codeLoc, requiredLoc) 
#runRates()


##### 
#  
#  Generate the lines of sight
#
#####
genLines(galID, gasfile, summaryLoc, expn, incline, nlos, maximpact, ncores)
print 'Lines generated'

##### 
#  
#  Run lines of sight through box with cellfinder
#
#####
command = './cellfinder'
print 'Cellfinder Done'


