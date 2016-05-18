
'''
Functions to genearte basic plots

Plots made:
    EW vs. impact parameter
    Ew distribution
    Covering fraction
    Column Density distribution
    Phase
'''

import covering as cv
import ewVsd as ew
import ewdist as ed
import nT as nt
import columndist as co

def make_plots(ions):

    # Make EW vs. Impact paramter       
    try:
        ew.ew_profile(ions)
    except:
        print 'Error running ewvsd'

    # Make EW distribution
    try:
        ed.ewdist(ions)
    except:
        print 'Error running ewdist'
    
    # Make covering fraction
    try:
        cv.covering_fraction(ions)
    except:
        print 'Error running covering'

    # Make column density distribuion
    try:
        co.column_distribution(ions)
    except:
        print 'Error running coldense'

    # Make phase plots
    try:
        nt.phase(ions)
    except:
        print 'Error running nT'

