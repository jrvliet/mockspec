
'''
Functions to genearte basic plots

Plots made:
    EW vs. impact parameter
    Ew distribution
    Covering fraction
    Column Density distribution
    Phase
'''

from . import covering as cv
from . import ewVsd as ew
from . import ewdist as ed
from . import nT as nt
from . import columndist as co

def make_plots(ions):

    # Make EW vs. Impact paramter       
    print('\tPlotting EW vs. Impact')
    try:
        ew.ew_profile(ions)
    except:
        print('Error running ewvsd')

    # Make EW distribution
    print('\tPlotting EW distribution')
    try:
        ed.ewdist(ions)
    except:
        print('Error running ewdist')
    ed.ewdist(ions)
    
    # Make covering fraction
    print('\tPlotting covering fraction')
    try:
        cv.covering_fraction(ions)
    except:
        print('Error running covering')

    # Make column density distribuion
    print('\tPlotting column density distribution')
    try:
        co.column_distribution(ions)
    except:
        print('Error running coldense')

    # Make phase plots
    print('\tPlotting phase')
    try:
        nt.phase(ions)
    except:
        print('Error running nT')

