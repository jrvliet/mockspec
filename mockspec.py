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



#  Read in the control file
props, ions, xh, instruments = readControlFile()
galID, expn, nlos, incline, ewcut, snr, ncores, rootLoc, requriedLoc = props

# Genearte the name of the gasfile
gasfile = galID+'_GZa'+expn+'.txt'
 

# Generate gal_props.dat file
setupGalprops( galID, expn )


##### 
#  
#  Run rates
#
#####
setupRatesControl( gasfile, expn, ions, requiredLoc)

rates(genLOS(galID, gasfile, summaryLoc, expn, incline, nlos, ncores)



##### 
#  
#  Generate the lines of sight
#
#####
genLOS(galID, gasfile, summaryLoc, expn, incline, nlos, ncores)


##### 
#  
#  Run lines of sight through box with cellfinder
#
#####
command = './cellfinder'



