
'''
Python version of cellfinder.c
Utilizes parallization, box cut, and dataframes
'''

from __future__ import print_function
import numpy as np
import pandas as pd
import joblib as jl

class los(object):

    '''
    Class to describe a line of sight vector
    '''
    
    def __init__ (self):

        # Entrance and exit points
        self.xen = 0
        self.yen = 0
        self.zen = 0
        self.xex = 0
        self.yex = 0
        self.zex = 0

        # Observable properties
        self.impact = 0
        self.phi = 0
        
        # Number
        self.number = 0

        # Directional Cosines
        self.losx = 0
        self.losy = 0
        self.losz = 0

        # Head of unit vector along LOS
        self.x0 = 0
        self.y0 = 0
        self.z0 = 0
        
        
def lmn(line):

    '''
    Calculates the directional cosines for the LOS
    '''

    # Determine the directional vector
    dx = line.xex - line.xen
    dy = line.yex - line.yen
    dz = line.zex - line.zen

    # Determine the directional cosines
    dl = np.sqrt(dx**2 + dy**2 + dz**2)
    line.losx = dx/dl
    line.losy = dy/dl
    line.losz = dz/dl
    
    # Determine a unit vector along the LOS
    line.x0 = line.xen + line.losx
    line.y0 = line.yen + line.losy
    line.z0 = line.zen + line.losz
    
    # Return the finished line
    return line
    

def cellsearch(line, gas, outfile):

    '''
    Determines the cells that lie along the line of sight defined by line
    '''

    # Get the minimum distance between each 
    gas['dmin'] = gas.apply(dmin_calc,line=line,axis=1)
    
    # Select out cells with dmin < distance between cell center and corner
    alongInds = gas['dmin']<gas['cell_size'].apply(lambda x: x*np.sqrt(3.)/2.)

    # Confirm the LOS actually passes through these cells
    gas['along'] = gas[alongInds].apply(confirm_cell,line=line,axis=1)
    gas.fillna(False,inplace=True,axis=1)
    #gas.to_csv('gasAlong.txt',sep=' ',index=False)
    
    cellsFound = gas['cellID'][gas['along']]

    # Write to file
    cellsFound.to_csv(outfile,index=False)


def confirm_cell(cell,line):

    '''
    Confirms the LOS actually passes through the cells
    '''

    # Equations of each wall
    xLeft = cell['x'] - cell['cell_size']/2.
    xRight = cell['x'] + cell['cell_size']/2.
    yLeft = cell['y'] - cell['cell_size']/2.
    yRight = cell['y'] + cell['cell_size']/2.
    zLeft = cell['z'] - cell['cell_size']/2.
    zRight = cell['z'] + cell['cell_size']/2.


    # Calculate intersection point
    d_x1 =(cell['x']-line.xen)**2 + (cell['y']-line.yen)**2 + (cell['z']-line.zen)**2    
    
    r_xpt = np.sqrt(d_x1 - cell['dmin']**2)
    
    # Find value of parameter 't' corresponding to the point r_xpt on the
    # parametric line
    #rdist = (cell['x']-line.x0)**2 + (cell['y']-line.y0)**2 + (cell['z']-line.z0)**2 
    #rXptObs = np.sqrt( rdist - cell['dmin']**2 )

    rxXpt = line.xen + line.losx*r_xpt
    ryXpt = line.yen + line.losy*r_xpt
    rzXpt = line.zen + line.losz*r_xpt

    # Does the point (rx_xpt,ry_xpt,rz_xpt) lie inside the cell?
    if rxXpt>xLeft and rxXpt<xRight and  ryXpt>yLeft and ryXpt<yRight and rzXpt>zLeft and rzXpt<zRight:
        return True
    else:
        return False
    

def dmin_calc(cell,line):

    ''' 
    Calculates the minimum distance between the LOS and the cell
    '''

    dx1 = cell['x'] - line.xen
    dy1 = cell['y'] - line.yen
    dz1 = cell['z'] - line.zen
    dx2 = cell['x'] - line.xex
    dy2 = cell['y'] - line.yex
    dz2 = cell['z'] - line.zex

    term1 = (dy1*dz2 - dz1*dy2)
    term2 = (dx1*dz2 - dz1*dx2)
    term3 = (dx1*dy2 - dy1*dx2)
    
    mag_xprod = np.sqrt(term1**2 + term2**2 + term3**2)
    mag_x2 = np.sqrt( (line.xex-line.xen)**2 +
                      (line.yex-line.yen)**2 +
                      (line.zex-line.zen)**2 )
    return mag_xprod/mag_x2


def cellfinder(run):


    # Read in gas box
    gasbox = '../{0:s}_GZa{1:s}.h5'.format(run.galID,run.expn)
    gas = pd.read_hdf(gasbox,'data')
    gas['cellID'] = [i+1 for i in gas.index.tolist()]
    
    # Convert the cellsize from parsecs to kiloparsecs
    gas['cell_size'] = gas['cell_size']/1000.
    maxCellSize = gas['cell_size'].max()

    # Read in LSO info
    linesHeader = 'xen yen zen xex yex zex'.split()
    losDat = pd.read_csv('lines.dat',skiprows=2,sep='\s+',names=linesHeader)
    
    linesHeader = 'losnum impact phi incline'.split()
    losInfo = pd.read_csv('lines.info',skiprows=2,sep='\s+',names=linesHeader)

    
    # Loop over lines of sight
    for i in range(run.nlos):

        # Generate los object
        line = los()
        line.xen = losDat['xen'].iloc[i]
        line.yen = losDat['yen'].iloc[i]
        line.zen = losDat['zen'].iloc[i]
        line.xex = losDat['xex'].iloc[i]
        line.yex = losDat['yex'].iloc[i]
        line.zex = losDat['zex'].iloc[i]
        line.impact = losInfo['impact'].iloc[i]
        line.phi = losInfo['phi'].iloc[i]
        line.number = i+1

        # Get the directional cosines (losx, losy, losz) and unit vector
        # describing the los vector through the box (tail at entrance point, 
        # head at (x0,y0,z0)
        line = lmn(line)

        # Select out a subset of box constrained by the entrance and exit points
        xPosLim = max([line.xen,line.xex]) + maxCellSize
        xNegLim = min([line.xen,line.xex]) - maxCellSize
        yPosLim = max([line.yen,line.yex]) + maxCellSize
        yNegLim = min([line.yen,line.yex]) - maxCellSize
        zPosLim = max([line.zen,line.zex]) + maxCellSize
        zNegLim = min([line.zen,line.zex]) - maxCellSize
        
        subbox = gas[ (gas['x']<xPosLim) & (gas['x']>xNegLim) & 
                      (gas['y']<yPosLim) & (gas['y']>yNegLim) & 
                      (gas['z']<zPosLim) & (gas['z']>zNegLim) ]

        print( len(gas), len(subbox), maxCellSize)

        # Generate filename 
        outfile = 'los{0:04d}.cellID.dat'.format(line.number)
        
        # Find all cells along the LOS
        cellsearch(line, subbox, outfile)

        













