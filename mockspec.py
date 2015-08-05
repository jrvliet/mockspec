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

##### 
#  
#  Set one: Generate the lines of sight
#
#####



