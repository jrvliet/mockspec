#
#  mockspec.py
# 
#  Full pipeline for performing synthetic QSO observations
#  of ART simulations of galaxy CGM

import numpy as np
import tables as tb
from mockspec_funcs import *
from genLOS import *



#  Read in the control file
props, ions, xh, instruments = readControlFile()
galID, expn, nlos, incline, ewcut, snr, ncores, rootLoc = props

# Genearte the name of the gasfile
gasfile = galID+'_GZa'+expn+'.txt'


##### 
#  
#  Step one: Generate the lines of sight
#
#####
genLOS(galID, gasfile, summaryLoc, expn, incline, nlos, ncores)



##### 
#  
#  Step two: Run lines of sight through box with cellfinder
#
#####

