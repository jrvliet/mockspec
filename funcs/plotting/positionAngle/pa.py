
'''
Functions to plot properites related to position angle
'''

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc


def convert_to_pa(phi):

    '''
    Function to convert a list of angles to position angle.
    Position angle is defined as the smallest angle between
    the LOS and the galaxy's major axis
    '''

    pa = []
    for angle in phi:

        if angle<=90:
            pa.append(angle)
        elif angle>90 and angle<180:
            pa.append( 180.0 - angle )
        elif angle>=180 and angle<360:
            pa.append( angle - 180.0 )
        else:
            pa.append( 360.0 - angle )
    return angle



def plot_imp_pa_binned(pa, b, value, stat, valLabel):

    '''
    Function to plot the mean EW of absorption at a given
    impact parameter and position angle
    '''

    numbins = 10

    # Bin the data and take the mean
    h, xedges, yedges, binnumber = sc.stats.binned_statistic_2d(pa, b, value, statistic=stat)
    h = np.ma.masked_where(h==0,h)

    # Plot
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(xedges, yedges, h)
    
    ax.set_xlabel('Position Angle [deg]')
    ax.set_ylabel('Impact Parameter [kpc]')

    cbar = plt.colorbar(mesh, ax=ax, use_gridspec=True)
    cbarLabel = '{0:s} {1:s}'.format(stat, valLabel)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel(cbarLabel, rotation=270, fontsize=12)

    plt.tight_layout()
    s = '{0:s}_{1:s}_imp_pa_bins.png'.format(valLabel, stat)
    plt.savefig(s, bbox_inches='tight')










    
